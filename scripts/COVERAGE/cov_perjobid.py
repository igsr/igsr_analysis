from PyHive.Hive import *
import subprocess
import pdb
import argparse
import re

parser = argparse.ArgumentParser(description='Calculate average cov for a job id in a Hive database')

parser.add_argument('--id', required=True, help='Job id to analyse' )

args = parser.parse_args()

db=HiveDB(host = 'mysql-rs-1kg-prod.ebi.ac.uk',
          user = 'g1krw',
          port = 4175,
          pwd = 'thousandgenomes',
          db = 'elowy_vc_1000g_highcov_30042019')

def print_coverage(jid):
    j=db.fetch_job_byid(id=int(jid))

    #get transpose bam url
    input_id=db.fetch_job_byid(j.prev_job_id).input_id
    pjobid=j.prev_job_id
    
    p = re.compile( '_extended_data_id \d+' )
    p1 = re.compile( '.*SQ_end.*' )
    m = p.match( input_id )
    input_id_reg=j.input_id
    while m is None:
        nj=db.fetch_job_byid(pjobid)
        input_id=nj.input_id
        if p1.match( input_id ): 
            input_id_reg=input_id
        pjobid=nj.prev_job_id
        m = p.match( input_id )

    (prefix,ad_id)=input_id.split(" ")
    ad=db.fetch_analysisdata_byid(id=int(ad_id))
    dataset,bam_url=ad.data.replace(" ","").split(",")[1].replace("\"","").replace("}","").split("=>")
    assert dataset=="bam", "No bam URL found"

    #get region
    region=input_id_reg.split(",")[-1].split("=>")[-1]
    region=region.replace("\"","").replace(" ","").replace("}","")
    (analysiset,chrom,start,chrom,end)=region.split(".")
    command = "samtools depth -a -d 0 -r %s:%s-%s %s | awk 'BEGIN {max = 0}"\
            "{if ($3>max) max=$3;sum+=$3;cnt++}END{print cnt \"\t\" sum \"\t\" max}'" \
            % (chrom, start, end, bam_url)
    p =  subprocess.Popen(command, stdout=subprocess.PIPE,shell=True)
    stdout,sterr = p.communicate()
    elms=stdout.decode("utf-8").split("\t")
    if elms[0]!="" and elms[1]!="": 
        bases_mapped, sum_of_depths, max= map(int, elms)
        avg=sum_of_depths/(int(end)-int(start))
        avg=round(avg, 2)
        print("{0} {1} {2} {3} {4} {5}\n".format(chrom,start,end,avg,j.job_id,j.status))
    else:
        print("{0} {1} {2} {3} {4} {5}\n".format(chrom,start,end,0,j.job_id,j.status))

print("#chr start end avg job_id status")
print_coverage(args.id)

