/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.survival.kaplanmeier.figure;

import org.biojava3.survival.data.WorkSheet;

/**
 *
 * @author willishf at gmail.com
 */
public interface CensorStatusSelect {

    /**
     *
     * @param worksheet
     * @param row
     * @return
     * @throws Exception
     */
    public CensorStatus select(WorkSheet worksheet,String row) throws Exception; 
}
