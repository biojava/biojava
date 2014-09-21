/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.phylo;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 *
 * @author willishf
 */
public class ProgressListenerStub implements NJTreeProgressListener {
	
	private static final Logger logger = LoggerFactory.getLogger(ProgressListenerStub.class);

	String priorState = "";
    public void progress(Object njtree,String state, int percentageComplete) {
        if(!priorState.equals(state)){
            priorState = state;
        }

        logger.info("{}:{}", state, percentageComplete);
    }



    public void progress(Object njtree,String state, int currentCount, int totalCount) {
        if (!priorState.equals(state)){
            priorState = state;
        }

        logger.info("{}:{}/", state, currentCount, totalCount);
    }

    public void complete(Object njtree) {

    }

    public void canceled(Object njtree) {
       
    }

}
