import os, glob
from utilities import cfg, job_runner
from utilities.job_runner import MyJob
from optparse import OptionParser

PACKAGE_DIR = "/".join(os.path.abspath(__file__).split("/")[:-2])

def parse_input_file(input_file):
    projs = {}
    
    for line in open(input_file, 'r'):
        cols = line.rstrip('\n').split()
        lib, proj = cols[0], cols[1]
        if not projs.has_key(proj):
            projs[proj] = []
        projs[proj].append(lib)
            
    return projs
        
def submit(cfg_file, input_file, options):
    task = 'polyA'
    script = '%s/CLEAT/polyA.py' % PACKAGE_DIR
    version = os.environ['TRANSABYSS_VERSION']
    
    jobs = []
    projs = parse_input_file(input_file)
    for proj, libs in projs.iteritems():
        for lib in libs:
            settings = cfg.initialize_settings(cfg_file, project=proj)
            
            top_dir = settings[proj]['topdir']
            working_dir = '%s/%s/trans-abyss-v%s' % (top_dir, lib, version)
            aln_dir = '%s/contigs_to_genome/%s-contigs/output' % (working_dir, lib)
            seq_dir = '%s/contigs_to_genome/%s-contigs/input' % (working_dir, lib)
            
            cluster_dir = '%s/CLEAT/cluster' % working_dir
            out_dir = '%s/CLEAT/tmp' % working_dir
            for d in (cluster_dir, out_dir):
                if not os.path.isdir(d):
                    os.makedirs(d)
            
            args = {}
            args['aln'] = '%s/seq.$TA_JOBID.sam' % aln_dir
            args['fa'] = '%s/seq.$TA_JOBID.fa' % seq_dir
            args['out'] = '%s/$TA_JOBID.tsv' % out_dir
            args['genome'] = settings[proj]['reference']
            args['bam'] = '%s/reads_to_contigs/%s-contigs.bam' % (working_dir, lib)
            
            if not os.path.isdir(seq_dir) or not os.path.isdir(aln_dir) or not os.path.exists(args['bam']):
                if not os.path.isdir(seq_dir):
                    print 'missing sequence dir:%s' % seq_dir
                if not os.path.isdir(aln_dir):
                    print 'missing sequence dir:%s' % aln_dir
                if not not os.path.exists(args['bam']):
                    print 'missing bam file:%s' % args['bam']
                sys.exit()
            
            seq_files = glob.glob(os.path.join(seq_dir, '*.fa'))
        
            job_name = '%s-polyA' % lib
            job_content = 'time python %s %s' % (script, cfg.create_args(settings, 'commands', task, args) )
            job = MyJob(content=job_content, name=job_name, config_key=task, first_task_id='1', last_task_id=str(len(seq_files)), queue=options.queue)
            job_runner.submit_jobs([job], settings, options.cluster, jobdir=cluster_dir, debug=False)
        
            submit_filter([job], lib, settings, working_dir, out_dir, cluster_dir, options)
            
def submit_filter(batch_job, lib, settings, working_dir, results_dir, cluster_dir, options):
    task = 'polyA-filter'
    script = '%s/CLEAT/polyA.py' % PACKAGE_DIR
    version = os.environ['TRANSABYSS_VERSION']
    
    args = {}
    args['results'] = results_dir
    args['out'] = '%s/CLEAT/polyA.tsv' % working_dir
    args['lib'] = lib
    
    job_name = '%s-polyA-filter' % lib
    job_content = 'time python %s %s' % (script, cfg.create_args(settings, 'commands', task, args) )
    job = MyJob(content=job_content, name=job_name, config_key=task, predecessors=batch_job, queue=options.queue)
    job_runner.submit_jobs([job], settings, options.cluster, jobdir=cluster_dir, debug=False)
    
    
def main(args, options):
    cfg_file, input_file = args
        
    submit(cfg_file, input_file, options)
    
if __name__ == '__main__':
    usage = "Usage: %prog config_file input_file"

    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--cluster", dest="cluster", help="cluster head node name", default='genesis')
    parser.add_option("-q", "--queue", dest="queue", help="cluster queue", default='all.q')
    
    (options, args) = parser.parse_args()
    main(args, options)
