$pdflatex = 'pdflatex --synctex=1 --shell-escape -interaction=nonstopmode';
sub build_header {
  system("ruby ./prep.rb")
}

build_header()