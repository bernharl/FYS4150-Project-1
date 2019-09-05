cd src

echo "Generate new data before plotting? (y/n)"
read yn
if [ "$yn" == "y" ]
then
  echo "Compiling c++ application used for main calculations"
  g++ -std=c++11  project1b.cpp -o project1b.out
  echo "Generating data using c++ application"
  ./project1b.out
fi


echo "Generating error plots"
python3 errorplot.py special.dat
echo "Generating cpu time plots"
python3 CPUtimeplot.py thomas.dat special.dat
echo "Compiling TeX report"
#cd ../doc/
#pdftex CompPhysProj1.tex
