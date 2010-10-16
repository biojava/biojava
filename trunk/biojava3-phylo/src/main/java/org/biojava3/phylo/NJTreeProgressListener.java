/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.phylo;

//import org.biojavax.bio.phylo.jalview.NJTree;


/**
 *
 * @author willishf
 */
public interface NJTreeProgressListener {
    public void progress(Object njtree,String state, int percentageComplete);
    public void progress(Object njtree,String state, int currentCount,int totalCount);
    public void complete(Object njtree);
    public void canceled(Object njtree);
}
