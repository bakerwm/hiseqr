---
title: "Demultiplexing report"
author: "Ming Wang"
date: "2020-06-16"
output:
  html_document:
    highlight: tango
    toc: yes
    toc_float:
      collapsed: no
    keep_md: true
  word_document:
    toc: yes
  pdf_document:
    toc: yes
params:
  input_dir: ""
---













## Summary


```
A total of 1 M reads in this lane, 0.8 M (0.8%) were assigned to 40 samples, 0.2 M (0.2%) were failed to assign.
```


## Parameters

+ Number of mismatches: 0 



## Table

<!--html_preserve--><div id="htmlwidget-d9bda3bdc38d4b5e1067" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d9bda3bdc38d4b5e1067">{"x":{"filter":"top","filterHTML":"<tr>\n  <td><\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"integer\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"2\" data-max=\"202905\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"0\" data-max=\"0.2\" data-scale=\"1\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41"],["ATACseq_Harwich_w1118_6h_rep1","ATACseq_Harwich_w1118_6h_rep2","ATACseq_Harwich_w1118_6h_rep3","ATACseq_Harwich_w1118_6h_rep4","ATACseq_w1118_Harwich_6h_rep1","ATACseq_w1118_Harwich_6h_rep2","ATACseq_w1118_Harwich_6h_rep3","ATACseq_w1118_Harwich_6h_rep4","DNAseq_Nulei_library_sgRNA_DOX_12_days_rep1","DNAseq_Nulei_library_sgRNA_DOX_12_days_rep2","DNAseq_Nulei_library_sgRNA_DOX_12_days_rep3","DNAseq_Nulei_library_sgRNA_DOX_12_days_rep4","DNAseq_Nulei_library_sgRNA_DOX_4_days_rep1","DNAseq_Nulei_library_sgRNA_DOX_4_days_rep2","DNAseq_Nulei_library_sgRNA_DOX_4_days_rep3","DNAseq_Nulei_library_sgRNA_DOX_4_days_rep4","DNAseq_Nulei_library_sgRNA_DOX_8_days_rep1","DNAseq_Nulei_library_sgRNA_DOX_8_days_rep2","DNAseq_Nulei_library_sgRNA_DOX_8_days_rep3","DNAseq_Nulei_library_sgRNA_DOX_8_days_rep4","DNAseq_Nulei_library_sgRNA_enrich_rep1","DNAseq_Nulei_library_sgRNA_enrich_rep2","DNAseq_Nulei_library_sgRNA_enrich_rep3","DNAseq_Nulei_library_sgRNA_enrich_rep4","DNAseq_pComb3XSS_nanobody","DNAseq_pEntry_nanobody","RNAseq_Harwich_w1118_6h_rep1","RNAseq_Harwich_w1118_6h_rep2","RNAseq_nosGal4XshCG4936_ovary_rep1","RNAseq_nosGal4XshCG4936_ovary_rep2","RNAseq_nosGal4XshWhite_ovary_rep1","RNAseq_nosGal4XshWhite_ovary_rep2","RNAseq_w1118_Harwich_6h_rep1","RNAseq_w1118_Harwich_6h_rep2","small_RNAseq_GFP_4_6h_IP_rep1","small_RNAseq_GFP_4_6h_IP_rep2","small_RNAseq_Piwi_4_6h_IP_rep1","small_RNAseq_Piwi_4_6h_IP_rep2","small_RNAseq_embryo_4_6h_input_rep1","small_RNAseq_embryo_4_6h_input_rep2","undemx"],[38198,36956,49323,27717,77252,54973,56410,70669,26,120,30594,30546,3985,10834,6081,30,2,12695,10,86,24129,44303,19689,18318,610,657,29666,23749,14837,22506,15684,33190,19734,23366,63,38,7,6,24,12,202905],[0,0,0,0,0.1,0.1,0.1,0.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>filename<\/th>\n      <th>count<\/th>\n      <th>million<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":10,"scrollX":true,"columnDefs":[{"className":"dt-right","targets":[2,3]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


**[#EOF]**
