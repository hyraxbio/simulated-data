pip uninstall simulate && pip install .
make-proviral-data --out /home/jeicher/Apps/simulated-data/_tests/output --sequences data/split/Seq7_Sus --platform roche --hypermutation-rate 3 --working-dir ./_tests/output/working
