/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * created at Sep 7, 2007
 */
package org.biojava.bio.structure.server;


import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import org.biojava.bio.structure.Structure;

/**
 * @deprecated
 * 
 *
 */
public class SimpleStructureServer 
implements StructureServer, 
StructureListener {

	public static final Logger logger = Logger.getLogger("org.biojava.bio.structure");

	List<Structure> queue;
	List<StructureListener> listeners; 
	PDBInstallation installation;

	int cacheSize;

	int countLoading;
	List<StructureListener> waitingList; 
	public  SimpleStructureServer(){

		listeners = new ArrayList<StructureListener>();
		queue = new ArrayList<Structure>();
		cacheSize = 1;
		countLoading = 0;
		waitingList = new ArrayList<StructureListener>();

	}

	public void initCache(){
		checkStatus();
	}

	public void addStructureListener(StructureListener listener) {
		listeners.add(listener);

	}

	public void clearStructureListeners() {
		listeners.clear();

	}

	public int getNrCPUs() {
		// TODO Auto-generated method stub
		return 0;
	}

	public PDBInstallation getPDBInstallation() {

		return installation; 
	}

	public synchronized void requestNextStructure(StructureListener listener) {

		if (installation ==null) {
			logger.warning("no PDB installation has been set, yet. Can not load next structure");
			triggerNextStructure(null, listener);
			return;
		}

		if (queue.size() > 0){

            // give the requestor a structure from the queue.

            Structure s = getStructureFromQueue();
			checkStatus();                        

            // and send the notification to requestor
			triggerNextStructure(s, listener);
			
		}

		// queue.size() == 0
		// on a multi CPU system this should not happen.
        // most likely we are on a single CPU system...
        // other possibility is that we are still loading the structures?

		if ( countLoading > 0 ){
			// will be the first one to be informed.
			
			waitingList.add(listener);
			
			checkStatus();
			return ;
		}

		System.err.println("queue is empty on server!");
		addNextToQueue(listener);
		checkStatus();
		return;

	}

    private synchronized void checkStatus(){
        System.out.println("server check status "+ queue.size() + "+" + countLoading + " == " + cacheSize + "?");
        while ( queue.size() + countLoading <= cacheSize ) {
            // refill the cache
            addNextToQueue(this);
            countLoading++;
        }

    }
    
	private void addNextToQueue(StructureListener listener){
		//System.out.println("server: addNextToQueue: queue:" + queue.size() + " " + installation.hasNext());
		if ( installation.hasNext()){

			StructureFetcherRunnable r = new StructureFetcherRunnable(installation);
			r.addStructureListener(listener);
			Thread t = new Thread(r);
			t.start();
		}
	}

	/** a new structure is loaded from the queue.
	 * remove it from queue and trigger loading of a new Structure 
	 * 
	 * @return Structure
	 */
	private synchronized Structure getStructureFromQueue(){
		if ( queue.size() == 0){
			return null;
		}
		Structure s = queue.get(0);
		queue.remove(0);		
		return s;
	}

	public void setCacheSize(int cacheSize) {
		this.cacheSize = cacheSize;

	}
	public int getCacheSize(){
		return cacheSize;
	}

	public void setPDBInstallation(PDBInstallation installation) {
		this.installation = installation;

	}

	public boolean hasNextStructure() {
		//System.out.println("server hasNextStructure loading:" +countLoading);
		if (queue.size() > 0)
			return true;
		if ( countLoading > 0)
			return true;
		boolean hasNext = installation.hasNext();
		if ( hasNext)
			return true;


		return false;
	}

	private void triggerNextStructure(Structure s, StructureListener li){
		StructureEvent event = new StructureEventImpl(s);

		li.newStructure(event);
	}



	public void modifiedStructure(StructureEvent event) {
		// TODO Auto-generated method stub

	}

	public synchronized void newStructure(StructureEvent event) {
		countLoading--;
		Structure s = event.getStructure();
		if ( s == null) {
			System.err.println("StructureServer: could not load structure for " + event.getPDBCode());
		} else {
			//System.out.println("server got new structure " + s.getPDBCode());
			if ( waitingList.size() > 0){
				StructureListener listener = waitingList.get(0);
                waitingList.remove(0);
				triggerNextStructure(s, listener);
				return;
			}
			// nobody is waiting for this structure, cache is
			queue.add(s);
		}

	}

	public void obsoleteStructure(StructureEvent event) {
		// TODO Auto-generated method stub

	}

}
