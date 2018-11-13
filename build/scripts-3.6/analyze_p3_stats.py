from p3BAMQC import p3BAMQC

print "h"

p=p3BAMQC(filepath="/Users/ernesto/projects/IGSR/files/20130606_sample_info.xlsx")

fqc_lc=p.get_final_qc_results(group="low coverage")
fqc_ex=p.get_final_qc_results(group="exome")

#getting samples failed chk_indel
failed_INDEL_lc_df=fqc_lc.loc[fqc_lc['Indel_Ratio']>5]
failed_INDEL_ex_df=fqc_ex.loc[fqc_ex['Indel_Ratio']>5]

header="#sample\tpopulation\tindel_ratio\tQC_passed\n"

f1=open('/Users/ernesto/projects/IGSR/files/failed_chkindel_lc.txt','w')
f1.write(header)
for index,row in failed_INDEL_lc_df.iterrows():
    newline="%s\t%s\t%s\t%s\n" % (index[0],index[1],row.Indel_Ratio,row.Passed_QC)
    f1.write(newline)
f1.close()

f2=open('/Users/ernesto/projects/IGSR/files/failed_chkindel_ex.txt','w')
f2.write(header)
for index,row in failed_INDEL_ex_df.iterrows():
    newline="%s\t%s\t%s\t%s\n" % (index[0],index[1],row.Indel_Ratio,row.Passed_QC)
    f2.write(newline)
f2.close()


####getting samples failing verifybamid
header1="#sample\tpopulation\tVerifyBam_Omni_Free\tVerifyBam_Affy_Free\tVerifyBam_Omni_Chip\tVerifyBam_Affy_Chip\tPassed_QC\n"

###LC
failed_VBAMID_lc_df=fqc_lc.loc[((fqc_lc.VerifyBam_Omni_Free>0.03) & (fqc_lc.VerifyBam_Omni_Chip>0.02)) | (fqc_lc.VerifyBam_Affy_Free>0.03) & ((fqc_lc.VerifyBam_Affy_Chip>0.02))]
f3=open('/Users/ernesto/projects/IGSR/files/failed_vbamid_lc.txt','w')
f3.write(header1)
for index,row in failed_VBAMID_lc_df.iterrows():
    newline="%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (index[0],index[1],row.VerifyBam_Omni_Free,row.VerifyBam_Affy_Free,row.VerifyBam_Omni_Chip,row.VerifyBam_Affy_Chip,row.Passed_QC)
    f3.write(newline)
f3.close()

all_failedQC_lc=fqc_lc.loc[fqc_lc.Passed_QC.isnull()]
f4=open('/Users/ernesto/projects/IGSR/files/all_failedQC_lc.txt','w')
f4.write(header1)
for index,row in all_failedQC_lc.iterrows():
    newline="%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (index[0],index[1],row.VerifyBam_Omni_Free,row.VerifyBam_Affy_Free,row.VerifyBam_Omni_Chip,row.VerifyBam_Affy_Chip,row.Passed_QC)
    f4.write(newline)
f4.close()


failed_VBAMID_ex_df=fqc_ex.loc[((fqc_ex.VerifyBam_Omni_Free>0.035) & (fqc_ex.VerifyBam_Omni_Chip>0.02)) | (fqc_ex.VerifyBam_Affy_Free>0.035) & ((fqc_ex.VerifyBam_Affy_Chip>0.02))]
f5=open('/Users/ernesto/projects/IGSR/files/failed_vbamid_ex.txt','w')
f5.write(header1)
for index,row in failed_VBAMID_ex_df.iterrows():
    newline="%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (index[0],index[1],row.VerifyBam_Omni_Free,row.VerifyBam_Affy_Free,row.VerifyBam_Omni_Chip,row.VerifyBam_Affy_Chip,row.Passed_QC)
    f5.write(newline)
f5.close()

all_failedQC_ex=fqc_ex.loc[fqc_ex.Passed_QC.isnull()]
f6=open('/Users/ernesto/projects/IGSR/files/all_failedQC_ex.txt','w')
f6.write(header1)
for index,row in all_failedQC_ex.iterrows():
    newline="%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (index[0],index[1],row.VerifyBam_Omni_Free,row.VerifyBam_Affy_Free,row.VerifyBam_Omni_Chip,row.VerifyBam_Affy_Chip,row.Passed_QC)
    f6.write(newline)
f6.close()

