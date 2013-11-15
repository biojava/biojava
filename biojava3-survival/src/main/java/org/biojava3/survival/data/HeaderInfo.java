/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.data;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class HeaderInfo {

    private Integer index;
    private boolean hide = false;

    /**
     *
     * @param index
     */
    public HeaderInfo(Integer index) {
        this.index = index;
    }

    @Override
    public String toString() {
        return index + ":" + hide;
    }

    /**
     * @return the index
     */
    public Integer getIndex() {
        return index;
    }

    /**
     * @param index the index to set
     */
    public void setIndex(Integer index) {
        this.index = index;
    }

    /**
     * @return the hide
     */
    public boolean isHide() {
        return hide;
    }

    /**
     * @param hide the hide to set
     */
    public void setHide(boolean hide) {
        this.hide = hide;
    }
    
    
}
