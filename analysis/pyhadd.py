import os,sys

addlist = ""
for n in xrange(400,500):
    f = "/net/nudsk0001/d00/scratch/taritree/scrape_out/output_scrape_%d.root"%(n)
    if not os.path.exists(f):
        continue
    stats = os.stat( f )
    print f," ",stats.st_size
    if stats.st_size>500:
        addlist += " %s"%(f)

os.system("hadd run400_499_scraped_allopt_tf3.root %s"%(addlist))
