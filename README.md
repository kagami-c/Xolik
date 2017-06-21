# Xolik
Finding cross-linked peptides with maximum paired scores in linear time

## Usage
Please ensure that the following packages exist.
- Python 2.7
- pandas

To run Xolik, use the following commands.
```
./Xolik -d db-filename.fasta -s msdata-filename.mzXML -o output-filename.csv --ms1tol 5 --ms2tol 0.02
python splitctrl.py output-filename.csv
```

Use `./Xolik --help` for all available options.

The "Score" field in the result files equals to -log10(e-value) by default,
and equals to XCorr if `--noevalue` is set.

## Dependencies
- TCLAP
- GoogleTest (optional)

## Publications
J. Dai\*, W. Jiang\*, F. Yu\* and W. Yu,
**"Xolik: finding cross-linked peptides with maximum paired scores in linear time"**,
in preparation. *Contributed equally to this work.

## License
BSD License