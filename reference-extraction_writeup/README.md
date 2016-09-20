See reference_extraction.pdf

Written in Pandoc Markdown, compiled with:

```
pandoc -S -s -f markdown -t latex --filter pandoc-fignos --filter=pandoc-citeproc -o $file_base_name.pdf -V mainfont:Garamond --latex-engine xelatex $file_name && google-chrome $file_base_name.pdf
```