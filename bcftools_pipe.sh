bcftools norm -m - ${1} | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bcftools view -e 'QUAL<30 || 0<ABHet<0.25' | bgzip -c > ${2}
