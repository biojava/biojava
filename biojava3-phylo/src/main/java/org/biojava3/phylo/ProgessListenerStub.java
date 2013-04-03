/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.phylo;


/**
 *
 * @author willishf
 */
public class ProgessListenerStub implements NJTreeProgressListener {
String priorState = "";
    public void progress(Object njtree,String state, int percentageComplete) {
        if(!priorState.equals(state)){
            priorState = state;
            System.out.println();
        }

        System.out.println("\n" + state + ":" + percentageComplete);
    }



    public void progress(Object njtree,String state, int currentCount, int totalCount) {
        if (!priorState.equals(state)){
            priorState = state;
            System.out.println();
        }

        System.out.println("\n" + state + ":" + currentCount + "/" + totalCount);
    }

    public void complete(Object njtree) {

    }

    public void canceled(Object njtree) {
       
    }

}
