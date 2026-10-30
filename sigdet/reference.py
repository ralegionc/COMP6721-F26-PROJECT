"""Our minimal reference-genome reader (no external dependencies).

To place a mutation in its 96-context bin we need the base just before and
just after it. Those come from the reference genome FASTA file, which is
about 3 GB, so we read it by random access using a .fai index: a small text
file listing where each chromosome starts in the FASTA. We build the index
once if it is missing (samtools faidx builds the same file faster).

The FASTA must be uncompressed (.fa/.fasta), and its build must match the
MAFs: GRCh38 for GDC files, GRCh37 for the MC3 public MAF.
"""
from collections import namedtuple
from pathlib import Path

# One .fai row: chromosome length, byte offset of its first base,
# bases per line, and bytes per line (bases plus the line ending).
FaiRecord = namedtuple("FaiRecord", "length offset line_bases line_bytes")


def build_fai(fasta_path, fai_path=None):
    """We scan the FASTA once and write a samtools-compatible .fai index."""
    fasta_path = Path(fasta_path)
    fai_path = Path(fai_path or str(fasta_path) + ".fai")
    rows = []
    # name/length/offset/lb/lw describe the current chromosome; pos is our byte position.
    name, length, offset, lb, lw, pos = None, 0, 0, None, None, 0
    with open(fasta_path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                # A header line starts a new chromosome; we save the previous one.
                if name is not None:
                    rows.append((name, length, offset, lb or 0, lw or 0))
                name = line[1:].split()[0].decode()
                length, lb, lw = 0, None, None
                offset = pos + len(line)
            else:
                seq = line.rstrip(b"\r\n")
                # The first sequence line fixes the line width for the whole chromosome.
                if lb is None and seq:
                    lb, lw = len(seq), len(line)
                length += len(seq)
            pos += len(line)
    if name is not None:
        rows.append((name, length, offset, lb or 0, lw or 0))
    fai_path.write_text("".join(f"{n}\t{l}\t{o}\t{b}\t{w}\n" for n, l, o, b, w in rows))
    return fai_path


class FastaReference:
    """We open the genome once and fetch short stretches of it by position."""

    def __init__(self, fasta_path):
        self.path = Path(fasta_path)
        if self.path.suffix == ".gz":
            raise ValueError("Decompress the reference first (gunzip); random access needs a plain FASTA")
        fai = Path(str(self.path) + ".fai")
        if not fai.exists():
            build_fai(self.path, fai)
        self.index = {}
        for line in fai.read_text().splitlines():
            n, l, o, b, w = line.split("\t")[:5]
            self.index[n] = FaiRecord(int(l), int(o), int(b), int(w))
        self._fh = open(self.path, "rb")

    def close(self):
        self._fh.close()

    # Lets us write "with FastaReference(path) as ref:" so the file always closes.
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def resolve(self, chrom):
        """We match MAF chromosome names to FASTA names ("chr1" vs "1", chrM vs MT)."""
        c = str(chrom)
        bare = c[3:] if c.lower().startswith("chr") else c
        candidates = [c, "chr" + bare, bare]
        if bare.upper() in {"M", "MT"}:
            candidates += ["chrM", "MT", "M"]
        for cand in candidates:
            if cand in self.index:
                return cand
        return None

    def fetch(self, chrom, start, end):
        """We return bases start..end (1-based, inclusive), uppercase. Positions
        outside the chromosome come back as N. An unknown chromosome returns ""."""
        name = self.resolve(chrom)
        if name is None:
            return ""
        rec = self.index[name]
        # Convert to 0-based, half-open coordinates clipped to the chromosome.
        s, e = max(start - 1, 0), min(end, rec.length)
        seq = ""
        if e > s:
            # Byte position of base i = offset + full lines before it + position within its line.
            a = rec.offset + (s // rec.line_bases) * rec.line_bytes + s % rec.line_bases
            last = e - 1
            b = rec.offset + (last // rec.line_bases) * rec.line_bytes + last % rec.line_bases + 1
            self._fh.seek(a)
            # The byte range may span a line break, which we strip out.
            seq = self._fh.read(b - a).replace(b"\n", b"").replace(b"\r", b"").decode().upper()
        left_pad = "N" * max(0, 1 - start)
        right_pad = "N" * max(0, end - rec.length)
        return left_pad + seq + right_pad
