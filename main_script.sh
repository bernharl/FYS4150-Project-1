cd src

echo "Generate new data before plotting? (y/n)"
read yn
if [ "$yn" == "y" ]
then
  echo "Compiling c++ application used for main calculations"
  g++ -std=c++11 -O2 project1b.cpp -o project1b.out
  echo "Generating data using c++ application"
  ./project1b.out
fi

echo "Generating function plot"
python3 plot1b.py data_thomas.dat data_special.dat

echo "Generating error plot"
python3 errorplot.py special.dat

echo "Generating cpu time plot"
python3 CPUtimeplot.py thomas.dat special.dat


echo "Build report? (y/n)"
read yn2
if [ "$yn2" == "y" ]
then
  cd ../doc/
  pdflatex -synctex=1 -interaction=nonstopmode CompPhysProj1.tex
fi
