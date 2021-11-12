python simulador.py -pE power.txt -pP coordinate.txt -n 200 -g 5 -s 1 -W 280 -H 280 -v 0 >log
python clusterizacaoEspectral.py -n 200 -g 5 -W 280 -H 280 -p power.txt -c coordinate.txt
python simulador.py -pE p.txt -pP c.txt -n 200 -G groupFile.txt -g 5 -s 1 -W 280 -H 280 -v 0 >log200nKME280m.txt
python SimParser.py -f log200nKME280m.txt >\python27\espectral-200nos-5grupos-280metros_s1.txt