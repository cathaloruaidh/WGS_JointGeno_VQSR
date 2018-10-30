#!/bin/bash

############################################################
#
# Assembly of NGS data from FASTQ to GVCF 
# Cathal Ormond 2017. 
#
# Alignment of reads from FASTQ file to reference genome
# Conversion to BAM file
#
# Based heavily on scripts found here:
# https://github.com/genepi-freiburg/gwas
#
############################################################




### Create a pdf report on the pipeline. 

if [[ -f ${RESULTS_DIR}/${REPORT_TEX} ]]
then
	rm ${RESULTS_DIR}/${REPORT_TEX}
fi


vt peek ${RESULTS_DIR}/${GENO_VCF_GZ} 2> ${RESULTS_DIR}/${GENO_VCF}.peek.txt
vt peek ${RESULTS_DIR}/${RECAL_VCF_GZ} 2> ${RESULTS_DIR}/${RECAL_VCF}.peek.txt
vt peek ${RESULTS_DIR}/${FILTER_VCF_GZ} 2> ${RESULTS_DIR}/${FILTER_VCF}.peek.txt


log "Annotation densities started." 3

${R} ${CURR_SCRIPT_DIR}/15_AnnotationDensities.R ${RESULTS_DIR}/${RECAL_TABLE} ${GRAPHICS_DIR}/${OUTPUT_FILE}


log "Tranches plots started. " 3

${R} ${CURR_SCRIPT_DIR}/15_TranchePlots.R ${RESULTS_DIR}/${SNP_TRANCHES} ${RESULTS_DIR}/${INDEL_TRANCHES} ${GRAPHICS_DIR}/${OUTPUT_FILE}

#exit 0


printf '\\documentclass[a4paper,11pt]{article}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\usepackage{amsmath, amsfonts, setspace, longtable, setspace, amssymb, fullpage, graphicx, pdfpages, underscore, fancyvrb}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\usepackage{grffile}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\usepackage[T1]{fontenc}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\usepackage[export]{adjustbox}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\input{/home/shared/tools/latex/macros.tex}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\parindent 0pt\n' >> ${RESULTS_DIR}/${REPORT_TEX}

printf '\\title{NGS Variant Recalibration Pipeline for %s}\n\n' $( echo ${OUTPUT_FILE} | sed -e 's/_/\\string_/g') >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\author{Cathal Ormond}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\date{%s}\n\n' "$(date '+%Y/%m/%d at %H:%M:%S')" >> ${RESULTS_DIR}/${REPORT_TEX}

printf '\\begin{document}\n\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\onehalfspacing\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\maketitle\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\tableofcontents\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\newpage\n\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}

printf '\\section{Post Joint-Genotyping Summary}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\begin{Verbatim}[fontsize=\scriptsize, frame=single, tabsize=4]\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '%s' "$(cat ${RESULTS_DIR}/${GENO_VCF}.peek.txt)" >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\end{Verbatim}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\newpage\n\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}

printf '\\section{Post Merge Summary}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\begin{Verbatim}[fontsize=\scriptsize, frame=single, tabsize=4]\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '%s' "$(cat ${RESULTS_DIR}/${RECAL_VCF}.peek.txt)" >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\end{Verbatim}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\newpage\n\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}

printf '\\section{Final Filtered VCF Summary}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\begin{Verbatim}[fontsize=\scriptsize, frame=single, tabsize=4]\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '%s' "$(cat ${RESULTS_DIR}/${FILTER_VCF}.peek.txt)" >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\end{Verbatim}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\newpage\n\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}

printf '\\section{\\texttt{peddy} Plots}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
for FILE in `ls ${GRAPHICS_DIR}/${OUTPUT_FILE}.recal.*png`
do
	MIDDLE=${FILE#${GRAPHICS_DIR}/${OUTPUT_FILE}.recal.}
	printf '\\subsection{%s}\n' $( echo ${MIDDLE%.png} | sed -e 's/_/\\string_/g' ) >> ${RESULTS_DIR}/${REPORT_TEX}
	printf '\\begin{center}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
	printf '\t\\includegraphics[max height=0.95\\textheight,max width=0.95\\textwidth]{%s}\n' $( echo ${FILE} | sed -e 's/_/\\string_/g' ) >> ${RESULTS_DIR}/${REPORT_TEX}
	printf '\\end{center}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
	printf '\\newpage\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
done

printf '\\section{VQSR Annotation Density Plots}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
for FILE in `ls ${GRAPHICS_DIR}/${OUTPUT_FILE}.recal.table.*pdf`
do
	MIDDLE=${FILE#${GRAPHICS_DIR}/${OUTPUT_FILE}.recal.table.}
	printf '\\subsection{%s}\n' ${MIDDLE%.pdf} >> ${RESULTS_DIR}/${REPORT_TEX}
	printf '\\begin{center}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
	printf '\t\\includegraphics[max height=0.95\\textheight,max width=0.95\\textwidth]{%s}\n' $( echo ${FILE} | sed -e 's/_/\\string_/g' ) >> ${RESULTS_DIR}/${REPORT_TEX}
	printf '\\end{center}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
	printf '\\newpage\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
done

printf '\\section{VQSR Tranches Plots}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\subsection{SNP Tranches Histogram \@ %s}\n' ${SNP_TRANCHE_LEVEL} >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\begin{center}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\t\\includegraphics[max height=0.95\\textheight,max width=0.95\\textwidth]{%s}\n' $( echo "${GRAPHICS_DIR}/${OUTPUT_FILE}.tranches.SNP.hist.pdf" | sed -e 's/_/\\string_/g' ) >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\end{center}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\subsection{SNP Ti/TV Ratio}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\begin{center}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\t\\includegraphics[max height=0.95\\textheight,max width=0.95\\textwidth]{%s}\n' $( echo "${GRAPHICS_DIR}/${OUTPUT_FILE}.tranches.SNP.line.pdf" | sed -e 's/_/\\string_/g' ) >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\end{center}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\subsection{INDEL Tranches Histogram \@ %s}\n' ${INDEL_TRANCHE_LEVEL} >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\begin{center}\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\t\\includegraphics[max height=0.95\\textheight,max width=0.95\\textwidth]{%s}\n' $( echo "${GRAPHICS_DIR}/${OUTPUT_FILE}.tranches.INDEL.hist.pdf" | sed -e 's/_/\\string_/g' ) >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\end{center}\n\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\n' >> ${RESULTS_DIR}/${REPORT_TEX}
printf '\\end{document}\n' >> ${RESULTS_DIR}/${REPORT_TEX}

pdflatex -output-directory ${RESULTS_DIR} ${RESULTS_DIR}/${REPORT_TEX} 
pdflatex -output-directory ${RESULTS_DIR} ${RESULTS_DIR}/${REPORT_TEX} \
#> ${LOG_DIR}/${OUTPUT_FILE}.ALIGN.${1}.log 2>&1

log "LaTeX" 3
