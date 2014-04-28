/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.kaplanmeier.figure;

import java.util.ArrayList;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GroupResults {


    /**
     *
     */
    public ArrayList<CensorStatus> group1;
    /**
     *
     */
    public ArrayList<CensorStatus> group2;
    /**
     *
     */
    public boolean group1WorseOutcome = false;
    /**
     *
     */
    public String signatureName ="";
}
