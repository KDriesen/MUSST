project: MUSST
summary: MUSST -- A modern Fortran MUltiScale Simulation Tribology solver<br/><br/> ![MUSST_img](media/MUSST_long_petit.png)
author:  Noël Brunetière - Arthur Francisco
email:   arthur.francisco@univ-poitiers.fr
website: https://www.pprime.fr/?q=fr/recherche-scientifique/d3/mecanique-des-interfaces-lubrifiees
github: https://github.com/Arthur-Francisco
project_github: https://github.com/TRIBO-Pprime/MUSST
src_dir: ./src
         ../DIGIS/src
         ../SPLIN/src
         ../UTILS/src
         ../MSOLV/src
include: ./inc
         ../MSOLV/inc
media_dir:  ./img
favicon:    ./img/logo.png
output_dir: ../docs
exclude_dir: ./src/arch
             ../SPLIN/src/arch/
exclude: prg.f90
         hsl_common90.f90
         hsl_ddeps90.f90
         hsl_ma48d.f90
fpp_extensions: F,f90
docmark:        !
docmark_alt:    #
predocmark:     >
predocmark_alt: <
display: public
         protected
         private
source: true
graph: true
search: true
css: ./css/MUSST.css
sort: src
coloured_edges: true
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc(anchorlink=False)
               markdown.extensions.smarty(smart_quotes=False)
---

-----------------
{!README.md!}

