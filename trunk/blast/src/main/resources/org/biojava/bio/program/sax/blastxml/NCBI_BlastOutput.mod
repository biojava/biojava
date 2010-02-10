<!-- ============================================ -->
<!-- This section mapped from ASN.1 module NCBI-BlastOutput -->
 
<!-- ============================================ -->
<!-- Definition of BlastOutput -->
 
 
<!--
$Id: NCBI_BlastOutput.mod,v 1.1 2003-06-25 02:04:18 dhuen Exp $
 
$Revision: 1.1 $
**********************************************************************
 
  ASN.1 for simplified BLAST output in XML
  by James Ostell, Yuri Wolf and Sergey Shavirin, 2000
 
 
 $Log: not supported by cvs2svn $
 Revision 6.6  2001/05/03 17:52:16  shavirin
 Adopted for usage with mani-iterational XML definition.
 
 Revision 6.4  2000/11/08 20:07:20  shavirin
 Added new parameter align_len analogos to the number reported in
 the Traditional Blast Output.
 
 Revision 6.3  2000/10/23 20:24:03  shavirin
 Few parameters were changed to OPTIONAL to allow XML without results:
 failure or no hits found condition.
 
 Revision 6.2  2000/08/11 17:48:35  shavirin
 Small fix.
 
 Revision 6.1  2000/08/09 20:43:12  shavirin
 Initial revision.
 
 
**********************************************************************
 -->
<!ELEMENT BlastOutput ( 
               BlastOutput_program ,
               BlastOutput_version ,
               BlastOutput_reference ,
               BlastOutput_db ,
               BlastOutput_query-ID ,
               BlastOutput_query-def ,
               BlastOutput_query-len ,
               BlastOutput_query-seq? ,
               BlastOutput_param ,
               BlastOutput_iterations )>
 
 
<!-- 
 BLAST program: blastp, tblastx etc.
 -->
<!ELEMENT BlastOutput_program ( #PCDATA )>
 
<!-- 
 Program version 
 -->
<!ELEMENT BlastOutput_version ( #PCDATA )>
 
<!-- 
 Steven, David, Tom and others
 -->
<!ELEMENT BlastOutput_reference ( #PCDATA )>
 
<!-- 
 BLAST Database name
 -->
<!ELEMENT BlastOutput_db ( #PCDATA )>
 
<!-- 
 SeqId of query
 -->
<!ELEMENT BlastOutput_query-ID ( #PCDATA )>
 
<!-- 
 Definition line of query
 -->
<!ELEMENT BlastOutput_query-def ( #PCDATA )>
 
<!-- 
 length of query sequence
 -->
<!ELEMENT BlastOutput_query-len ( %INTEGER; )>
 
<!-- 
 query sequence itself
 -->
<!ELEMENT BlastOutput_query-seq ( #PCDATA )>
 
<!-- 
 search parameters
 -->
<!ELEMENT BlastOutput_param ( Parameters )>
<!ELEMENT BlastOutput_iterations ( Iteration+ )>
 
 
 
 
<!-- Definition of Iteration -->
 
<!ELEMENT Iteration ( 
               Iteration_iter-num ,
               Iteration_hits? ,
               Iteration_stat? ,
               Iteration_message? )>
 
 
<!-- 
 iteration number
 -->
<!ELEMENT Iteration_iter-num ( %INTEGER; )>
 
<!-- 
 Hits one for every db sequence
 -->
<!ELEMENT Iteration_hits ( Hit* )>
 
<!-- 
 search statistics            
 -->
<!ELEMENT Iteration_stat ( Statistics )>
 
<!-- 
 Some (error?) information
 -->
<!ELEMENT Iteration_message ( #PCDATA )>
 
 
 
 
<!-- Definition of Parameters -->
 
<!ELEMENT Parameters ( 
               Parameters_matrix? ,
               Parameters_expect ,
               Parameters_include? ,
               Parameters_sc-match? ,
               Parameters_sc-mismatch? ,
               Parameters_gap-open ,
               Parameters_gap-extend ,
               Parameters_filter? ,
               Parameters_pattern? ,
               Parameters_entrez-query? )>
 
 
<!-- 
 Matrix used (-M)
 -->
<!ELEMENT Parameters_matrix ( #PCDATA )>
 
<!-- 
 Expectation threshold (-e)
 -->
<!ELEMENT Parameters_expect ( %REAL; )>
 
<!-- 
 Inclusion threshold (-h)
 -->
<!ELEMENT Parameters_include ( %REAL; )>
 
<!-- 
 match score for NT (-r)
 -->
<!ELEMENT Parameters_sc-match ( %INTEGER; )>
 
<!-- 
 mismatch score for NT (-q)
 -->
<!ELEMENT Parameters_sc-mismatch ( %INTEGER; )>
 
<!-- 
 Gap opening cost (-G)
 -->
<!ELEMENT Parameters_gap-open ( %INTEGER; )>
 
<!-- 
 Gap extension cost (-E)
 -->
<!ELEMENT Parameters_gap-extend ( %INTEGER; )>
 
<!-- 
 Filtering options (-F)
 -->
<!ELEMENT Parameters_filter ( #PCDATA )>
 
<!-- 
 PHI-BLAST pattern
 -->
<!ELEMENT Parameters_pattern ( #PCDATA )>
 
<!-- 
 Limit of request to Entrez query
 -->
<!ELEMENT Parameters_entrez-query ( #PCDATA )>
 
 
<!-- Definition of Statistics -->
 
<!ELEMENT Statistics ( 
               Statistics_db-num ,
               Statistics_db-len ,
               Statistics_hsp-len ,
               Statistics_eff-space ,
               Statistics_kappa ,
               Statistics_lambda ,
               Statistics_entropy )>
 
 
<!-- 
 Number of sequences in BLAST db
 -->
<!ELEMENT Statistics_db-num ( %INTEGER; )>
 
<!-- 
 Length of BLAST db
 -->
<!ELEMENT Statistics_db-len ( %INTEGER; )>
 
<!-- 
 Effective HSP length
 -->
<!ELEMENT Statistics_hsp-len ( %INTEGER; )>
 
<!-- 
 Effective search space
 -->
<!ELEMENT Statistics_eff-space ( %REAL; )>
 
<!-- 
 Karlin-Altschul parameter K
 -->
<!ELEMENT Statistics_kappa ( %REAL; )>
 
<!-- 
 Karlin-Altschul parameter Lambda
 -->
<!ELEMENT Statistics_lambda ( %REAL; )>
 
<!-- 
 Karlin-Altschul parameter H
 -->
<!ELEMENT Statistics_entropy ( %REAL; )>
 
 
<!-- Definition of Hit -->
 
<!ELEMENT Hit ( 
               Hit_num ,
               Hit_id ,
               Hit_def ,
               Hit_accession ,
               Hit_len ,
               Hit_hsps? )>
 
 
<!-- 
 hit number
 -->
<!ELEMENT Hit_num ( %INTEGER; )>
 
<!-- 
 SeqId of subject
 -->
<!ELEMENT Hit_id ( #PCDATA )>
 
<!-- 
 definition line of subject
 -->
<!ELEMENT Hit_def ( #PCDATA )>
 
<!-- 
 accession
 -->
<!ELEMENT Hit_accession ( #PCDATA )>
 
<!-- 
 length of subject
 -->
<!ELEMENT Hit_len ( %INTEGER; )>
 
<!-- 
 all HSP regions for the given subject
 -->
<!ELEMENT Hit_hsps ( Hsp* )>
 
 
 
<!-- Definition of Hsp -->
 
<!ELEMENT Hsp ( 
               Hsp_num ,
               Hsp_bit-score ,
               Hsp_score ,
               Hsp_evalue ,
               Hsp_query-from ,
               Hsp_query-to ,
               Hsp_hit-from ,
               Hsp_hit-to ,
               Hsp_pattern-from? ,
               Hsp_pattern-to? ,
               Hsp_query-frame? ,
               Hsp_hit-frame? ,
               Hsp_identity? ,
               Hsp_positive? ,
               Hsp_gaps? ,
               Hsp_align-len? ,
               Hsp_density? ,
               Hsp_qseq ,
               Hsp_hseq ,
               Hsp_midline? )>
 
 
<!-- 
 HSP number
 -->
<!ELEMENT Hsp_num ( %INTEGER; )>
 
<!-- 
 score (in bits) of HSP
 -->
<!ELEMENT Hsp_bit-score ( %REAL; )>
 
<!-- 
 score of HSP
 -->
<!ELEMENT Hsp_score ( %REAL; )>
 
<!-- 
 e-value of HSP
 -->
<!ELEMENT Hsp_evalue ( %REAL; )>
 
<!-- 
 start of HSP in query
 -->
<!ELEMENT Hsp_query-from ( %INTEGER; )>
 
<!-- 
 end of HSP
 -->
<!ELEMENT Hsp_query-to ( %INTEGER; )>
 
<!-- 
 start of HSP in subject
 -->
<!ELEMENT Hsp_hit-from ( %INTEGER; )>
 
<!-- 
 end of HSP in subject
 -->
<!ELEMENT Hsp_hit-to ( %INTEGER; )>
 
<!-- 
 start of PHI-BLAST pattern
 -->
<!ELEMENT Hsp_pattern-from ( %INTEGER; )>
 
<!-- 
 end of PHI-BLAST pattern
 -->
<!ELEMENT Hsp_pattern-to ( %INTEGER; )>
 
<!-- 
 translation frame of query
 -->
<!ELEMENT Hsp_query-frame ( %INTEGER; )>
 
<!-- 
 translation frame of subject
 -->
<!ELEMENT Hsp_hit-frame ( %INTEGER; )>
 
<!-- 
 number of identities in HSP
 -->
<!ELEMENT Hsp_identity ( %INTEGER; )>
 
<!-- 
 number of positives in HSP
 -->
<!ELEMENT Hsp_positive ( %INTEGER; )>
 
<!-- 
 number of gaps in HSP
 -->
<!ELEMENT Hsp_gaps ( %INTEGER; )>
 
<!-- 
 length of the alignment used
 -->
<!ELEMENT Hsp_align-len ( %INTEGER; )>
 
<!-- 
 score density
 -->
<!ELEMENT Hsp_density ( %INTEGER; )>
 
<!-- 
 alignment string for the query (with gaps)
 -->
<!ELEMENT Hsp_qseq ( #PCDATA )>
 
<!-- 
 alignment string for subject (with gaps)
 -->
<!ELEMENT Hsp_hseq ( #PCDATA )>
 
<!-- 
 formating middle line
 -->
<!ELEMENT Hsp_midline ( #PCDATA )>
 
 
 
 
