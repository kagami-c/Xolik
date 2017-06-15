# Xolik
Finding cross-linked peptides with maximum paired scores in linear time

## Usage
```bash
./Xolik -d db-filename.fasta -s msdata-filename.mzXML -o output-filename.csv --ms1tol 5 --ms2tol 0.02 --parallel --thread 4
```

Please use `./Xolik --help` for all available options.

## Dependencies
- TCLAP
- GoogleTest (optional)

## Publications
J. Dai\*, W. Jiang\*, F. Yu\* and W. Yu,
**"Xolik: finding cross-linked peptides with maximum paired scores in linear time"**,
in preparation. *Joint first authors.

## License
BSD License