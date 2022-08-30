for i in $(find . -iname "*.cc"); do
  git mv "$i" "${i/.cc/.cpp}";
done
