class QueryObj(object):
    '''
    todofile=None
    max_downloads=5
    num_streams=1
    donefile=~/git/scqc/test/query-donelist.txt
    
    '''
    def __init__(self, config):
        self.log = logging.getLogger('query')
        self.log.info('query init...')
        self.config = config
        self.todofile = self.config.get('query','todofile')
        if self.todofile.lower().strip() == "none":
            self.todofile = None
        else:
            self.todofile = os.path.expanduser(self.todofile)         
        
        self.donefile = self.config.get('query','donefile')
        if self.donefile.lower().strip() == "none":
            self.donefile = None
        else:
            self.donefile = os.path.expanduser(self.donefile)    
        self.shutdown = False
        self.sleep = int(self.config.get('query','sleep'))
        self.outlist = [] 
        
    def run(self):
        self.log.info('query run...')
        try:
            while not self.shutdown:
                self.log.debug(f'query cycle. {self.sleep} seconds...')
                self.todolist = self.readlist(self.todofile)
                self.donelist = self.readlist(self.donefile)
                if self.todolist is not None:
                    self.dolist = listdiff(self.todolist, self.donelist) 
                else:
                    self.dolist = None
                
                self.finished = self.execute(self.dolist)
                
                self.log.debug(f"got donelist of length {len(self.donelist)}. writing...")
                #self.writedone(self.finished)
                self.adddone(self.finished)
                self.log.debug(f"done writing donelist: {self.donefile}. sleeping {self.sleep} ...")
                time.sleep(self.sleep)
        
        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')
            
        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.log.error(traceback.format_exc(None))
            raise ex        
    
    def execute(self, todolist):
        '''
        Perform one run. 
        
        '''
        self.log.debug('executing...')
        self.log.info('ignoring todolist for query.')       
        sq = sra.Query(self.config)
        outlist = sq.execute()
        return outlist
        
          
    def stop(self):
        self.log.info('stopping...')        

    def readlist(self, filepath):
        if filepath is not None:
            self.log.info(f'reading file: {filepath}')
            flist = []
            try:
                with open(self.todofile, 'r') as f:
                   flist = [line.strip() for line in  f]
                self.log.debug(f'got todolist with {len(flist)} items.')
                return list 
            except:
                pass
        else:
            self.log.info('no file. return None.')
            return None


    def writedone(self, donelist):
        if self.donefile is not None:
            self.log.info('writing donefile...')
            rootpath = os.path.dirname(self.donefile)
            basename = os.path.basename(self.donefile)
            try:
                (tfd, tfname) = tempfile.mkstemp(suffix=None, 
                                              prefix=f"{basename}.", 
                                              dir=f"{rootpath}/", 
                                              text=True)
                self.log.debug(f"made temp {tfname}")
                with os.fdopen(tfd, 'w') as f:
                    nlines = 0
                    for item in donelist:
                        f.write(f"{item}\n")
                        nlines += 1
                os.rename(tfname, self.donefile)
                self.log.info(f"wrote {nlines} to {self.donefile}")
            except Exception as ex:
                self.log.error(traceback.format_exc(None))
                
            finally:
                pass
                #tfd.close()
                #os.remove(tfname, dir_fd=None)
        else:
            self.log.info('no donefile defined.')


    def listdiff(self, list1, list2):
        s1 = set(list1)
        s2 = set(list2)
        sd = s1 - s2
        dl = list(sd)
        dl.sort()
        return dl

class DownloadObj(object):
    '''
    todofile=~/git/scqc/test/query-donelist.txt
    max_downloads=5
    num_streams=1
    donelist=~/git/scqc/test/download-donelist.txt
    
    '''
    def __init__(self, config):
        self.log = logging.getLogger('download')
        self.log.info('download init...')
        self.config = config
        self.todofile = self.config.get('download','todofile')
        if self.todofile.lower().strip() == "none":
            self.todofile = None
        else:
            self.todofile = os.path.expanduser(self.todofile)         
        
        self.donefile = self.config.get('download','donefile')
        if self.donefile.lower().strip() == "none":
            self.donefile = None
        else:
            self.donefile = os.path.expanduser(self.donefile)   
         
        self.shutdown = False
        self.sleep = int(self.config.get('download','sleep'))
        self.outlist = [] 
        self.max_downloads = int(self.config.get('download','max_downloads'))
        self.num_streams = int(self.config.get('download','num_streams'))


        
    def run(self):
        self.log.info('downloader run...')
        try:
            while not self.shutdown:
                self.log.debug(f'download cycle. {self.sleep} seconds...')
                self.todolist = self.readlist(self.todofile)
                self.donelist = self.readlist(self.donefile)
                self.dolist = listdiff(self.todolist, self.donelist)
                self.donelist = self.execute(self.dolist)               
                self.log.debug(f"got donelist of length {len(self.donelist)}. writing...")
                self.writedone(self.donelist)
                self.log.debug(f"done writing donelist: {self.donefile}. sleeping {self.sleep} ...")
                time.sleep(self.sleep)
        
        except KeyboardInterrupt:
            print('\nCtrl-C. stopping.')
            
        except Exception as ex:
            self.log.warning("exception raised during main loop.")
            self.log.error(traceback.format_exc(None))
            raise ex        

    def execute(self, todolist):
        '''
        Perform one run. 
        
        '''
        self.log.debug('executing...')
        
        
        
        
        outlist = []
        dq = Queue()
        for runid in todolist:
            pf = sra.Prefetch(self.config, runid, outlist)
            dq.put(pf)
        logging.debug(f'created queue of {dq.qsize()} items')
        md = int(self.config.get('sra','max_downloads'))
        for n in range(md):
            sra.Worker(dq).start()
        logging.debug('waiting to join threads...')
        dq.join()
        logging.debug('all workers done...')
        return outlist


    
    def stop(self):
        self.log.info('stopping...')        

    def getdone(self):
        '''
            Get definitive list of already-completed ids. 
        '''
        
        
        

    def readlist(self, filepath):
        if filepath is not None:
            self.log.info(f'reading file: {filepath}')
            flist = []
            try:
                with open(self.todofile, 'r') as f:
                   flist = [line.strip() for line in  f]
                self.log.debug(f'got todolist with {len(flist)} items.')
                return list 
            except:
                pass
        else:
            self.log.info('no file. return None.')
            return None


    def writedone(self, donelist):
        if self.donefile is not None:
            self.log.info('writing donefile...')
            rootpath = os.path.dirname(self.donefile)
            basename = os.path.basename(self.donefile)
            try:
                (tfd, tfname) = tempfile.mkstemp(suffix=None, 
                                              prefix=f"{basename}.", 
                                              dir=f"{rootpath}/", 
                                              text=True)
                self.log.debug(f"made temp {tfname}")
                with os.fdopen(tfd, 'w') as f:
                    nlines = 0
                    for item in donelist:
                        f.write(f"{item}\n")
                        nlines += 1
                os.rename(tfname, self.donefile)
                self.log.info(f"wrote {nlines} to {self.donefile}")
            except Exception as ex:
                self.log.error(traceback.format_exc(None))
                
            finally:
                pass
                #tfd.close()
                #os.remove(tfname, dir_fd=None)
        else:
            self.log.info('no donefile defined.')
        