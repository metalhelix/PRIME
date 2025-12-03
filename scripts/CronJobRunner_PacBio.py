import logging
import os
from subprocess import check_output, CalledProcessError, STDOUT, call
import argparse
import pandas as pd
from pathlib import Path

def system_call(command: list) -> tuple:
    try:
        output = check_output(command, shell = True, stderr=STDOUT).decode()
        success = True
    except CalledProcessError as e:
        output = e.output.decode()
        success = False
    return output, success

def check_log(xml_paths: list, logfile): # log the whole path of the metadata.xml file! 
    """
    Check log file if fc id is already there, don't add if it is
    """
    seen_key = []
    with open(logfile, 'r') as f:
        for line in f:
            for xml in xml_paths:
                if xml in line:
                    seen_key.append(xml)
    xml_dict = {val for val in xml_paths if val not in seen_key}
    print(xml_dict)
    return xml_dict

def new_fc_pipeline_run(logfile, data_folder):

    logging.basicConfig(filename=logfile,
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')

    output, success = system_call([
            f"/n/ngs/tools/PRIME/PRIME/bin/fd '*\.metadata\.xml' \
                {data_folder} \
                -d 4 --glob --changed-within 128h --type f |grep -v 'basic_preview' |grep -v 'full_preview'"
        ])
    if not success:
            print("Something isn't working! Check the code")

    xml_paths = [out for out in output.split("\n") if os.path.dirname(out) != '']
    xml_dict = check_log(xml_paths, logfile) # Remove fc from dict if it's alraedy in the log file
    
    if len(xml_dict) > 0:
        for xml in xml_dict:
            logging.info(f"Discovered {xml}")
            os.chdir(Path(xml).parent.parent) # Change directory to the flowcell folder
            fcid = str(Path(xml).parent.parent).split("/")[-2] + "-" + str(Path(xml).parent.parent).split("/")[-1]
            cmd = f"source ~/.bash_profile; source /etc/profile; source ~/.bashrc; PATH=\"/usr/local/bin:$PATH\"; \
                    ml nextflow; nextflow run /n/core/Bioinformatics/analysis/CompBio/boweny/nf-Pipeline/PRIME/main.nf \
                    --PB_xml {xml} > \
                    /n/ngs/tools/PRIME/logs/machine_run_logs/{fcid}.PACBIO.nf.log"
                    
            recipients = "bioinfo@stowers.org mpe@stowers.org hhassan@stowers.org mcm@stowers.org dw2733@stowers.org"
            subject = f"Flowcell: {fcid}; Order Type: PacBio; Machine: Revio; PRIME execution STARTED"

            # Create the email content as a heredoc
            email_body = f"""/usr/sbin/sendmail -t <<EOF
To: {recipients}
From: bioinfo@stowers.org
Subject: {subject}
Content-Type: text/plain

The PRIME pipeline has STARTED for flowcell {fcid}.
EOF
"""

            # Send the email
            call(email_body, shell=True)

            # Run the command
            print(cmd)
            return_code = call(cmd, shell=True)
            if return_code == 0:
                logging.info(f"PRIME Nextflow Pipeline executed successfully {fcid}")
                print("Command executed successfully")
            else:
                print("Command failed with exit code", return_code)
            
    else:
        logging.info("No new runs")

PACBIO_log = "/n/ngs/tools/PRIME/logs/PACBIO_Orders.log"
## Check PacBio runs ###
data_folder = "/n/analysis/technology/revio/"
new_fc_pipeline_run(PACBIO_log, data_folder) 
