src_dir: ./src
         ../DIGIS/src
         ../SPLIN/src
         ../UTILS/src
         ../MSOLV/src
include: ./inc
         ../MSOLV/inc
media_dir:  ./img
output_dir: ./doc
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
source: false
graph: false
search: false
css: ./css/MUSST.css
