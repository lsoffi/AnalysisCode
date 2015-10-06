# Take a list of root file in txt file as input
# and make a header file as output
# containing a vector of TString

# $1: list.txt
# $2: tag

echo "vector<TString> list_"$2"()"
echo "{"
echo "vector<TString> list;"

grep root $1 | awk '{print "list.push_back(\042dcap://maite.iihe.ac.be/"$1"\042);"}'

echo "return list;"
echo "}"
