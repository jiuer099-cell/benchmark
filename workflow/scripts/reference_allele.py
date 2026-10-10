from pathlib import Path

class Reference:
    def __init__(self, path):
        self.stream = Path(path).open("rb")
        with Path(str(path) + ".fai").open() as handle:
            self.index = {f[0]: tuple(map(int, f[1:5])) for f in (line.rstrip().split("\t") for line in handle)}

    def close(self):
        self.stream.close()

    def fetch(self, chrom, start, end):
        length, offset, bases, width = self.index[chrom]
        if not 0 <= start <= end <= length:
            raise ValueError("reference interval is out of bounds")
        if start == end:
            return ""
        first = offset + start // bases * width + start % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases + 1
        self.stream.seek(first)
        result = self.stream.read(last - first).replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
        if len(result) != end - start:
            raise ValueError("reference FASTA/index mismatch")
        return result

    def equivalent(self, a, b):
        # Do not infer equivalence from graph IDs alone, size, or proximity.
        if a[0] != b[0]:
            return False
        left = min(int(a[1]), int(b[1])) - 1
        right = max(int(a[1]) - 1 + len(a[3]), int(b[1]) - 1 + len(b[3]))
        if right - left > 1_000_000:
            return False
        seq = self.fetch(a[0], left, right)
        haplotypes = []
        for f in (a, b):
            offset = int(f[1]) - 1 - left
            if seq[offset:offset + len(f[3])] != f[3].upper():
                return False
            if not set(f[3].upper() + f[4].upper()) <= set("ACGT"):
                return False
            haplotypes.append(seq[:offset] + f[4].upper() + seq[offset + len(f[3]):])
        return haplotypes[0] == haplotypes[1]


