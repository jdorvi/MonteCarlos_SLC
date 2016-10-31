declare -a files=("44017h2002.txt" "44017h2003.txt" "44017h2004.txt" "44017h2005.txt" "44017h2006.txt" "44017h2007.txt" "44017h2008.txt" "44017h2009.txt" "44017h2010.txt" "44017h2011.txt" "44017h2013.txt" "44017h2014.txt" "44017h2015.txt")
for filed in "${files[@]}"
do
    python /media/sf_C_DRIVE/Users/jdorvinen/Documents/Jared/Projects/East\ Hampton/met_data/montauk/met_data_10312016.py $filed
done

