/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.survival.kaplanmeier.metadata;

import org.biojava3.survival.data.WorkSheet;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public interface DiscreteQuantizerInterface {
    /**
     *
     * @param worksheet
     * @param column
     */
    public void process(WorkSheet worksheet, String column);

}
