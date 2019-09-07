cd src
echo "Generate new data before plotting? (y/n)"
read yn
if [ "$yn" == "y" ] # If y, compile and run both c++ codes with O3 optimization
then
  echo "Compiling c++ application used for main calculations"
  g++ -std=c++11 -O3 project1b.cpp -o project1b.out
  echo "Generating data using compiled application"
  ./project1b.out
  echo "Compiling LU code"
  g++ -std=c++11 -O3 project1e.cpp -larmadillo -o project1e.out
  echo "Generating data using compiled LU application"
  ./project1e.out

fi

 # Run plot script with corresponding data files
echo "Generating function plot"
python3 plot1b.py data_special_%d.dat data_thomas_%d.dat

echo "Generating error plot"
python3 errorplot.py special.dat

echo "Generating cpu time plot"
python3 CPUtimeplot.py thomas.dat special.dat


echo "Build report? (y/n)"
read yn2
# If y, compile TeX document. The compilation is run many times because
# bibtex is usually non-cooperative...
if [ "$yn2" == "y" ]
then
  cd ../doc/
  pdflatex -synctex=1 -interaction=nonstopmode CompPhysProj1.tex
  bibtex CompPhysProj1.aux
  pdflatex -synctex=1 -interaction=nonstopmode CompPhysProj1.tex
  bibtex CompPhysProj1.aux
  pdflatex -synctex=1 -interaction=nonstopmode CompPhysProj1.tex
fi
