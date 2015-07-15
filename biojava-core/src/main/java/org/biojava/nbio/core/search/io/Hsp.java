/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.core.search.io;

/**
 *
 * @author pavanpa
 */
public abstract class Hsp {
    private final int hspNum;
    private final double hspBitScore;
    private final int hspScore;
    private final double hspEvalue;
    private final int hspQueryFrom;
    private final int hspQueryTo;
    private final int hspHitFrom;
    private final int hspHitTo;

    public int getHspNum() {
        return hspNum;
    }

    public double getHspBitScore() {
        return hspBitScore;
    }

    public int getHspScore() {
        return hspScore;
    }

    public double getHspEvalue() {
        return hspEvalue;
    }

    public int getHspQueryFrom() {
        return hspQueryFrom;
    }

    public int getHspQueryTo() {
        return hspQueryTo;
    }

    public int getHspHitFrom() {
        return hspHitFrom;
    }

    public int getHspHitTo() {
        return hspHitTo;
    }

    public int getHspQueryFrame() {
        return hspQueryFrame;
    }

    public int getHspHitFrame() {
        return hspHitFrame;
    }

    public int getHspIdentity() {
        return hspIdentity;
    }

    public int getHspPositive() {
        return hspPositive;
    }

    public int getHspGaps() {
        return hspGaps;
    }

    public int getHspAlignLen() {
        return hspAlignLen;
    }
    /**
     * HSP aligned query sequence string
     * @return 
     */
    public String getHspQseq() {
        return hspQseq;
    }
    /**
     * HSP aligned hit sequence string
     * @return 
     */
    public String getHspHseq() {
        return hspHseq;
    }
    /**
     * Identity string representing correspondence between aligned residues
     * @return 
     */
    public String getHspIdentityString() {
        return hspIdentityString;
    }
    
    

    public Hsp(int hspNum, double hspBitScore, int hspScore, double hspEvalue, int hspQueryFrom, int hspQueryTo, int hspHitFrom, int hspHitTo, int hspQueryFrame, int hspHitFrame, int hspIdentity, int hspPositive, int hspGaps, int hspAlignLen, String hspQseq, String hspHseq, String hspIdentityString) {
        this.hspNum = hspNum;
        this.hspBitScore = hspBitScore;
        this.hspScore = hspScore;
        this.hspEvalue = hspEvalue;
        this.hspQueryFrom = hspQueryFrom;
        this.hspQueryTo = hspQueryTo;
        this.hspHitFrom = hspHitFrom;
        this.hspHitTo = hspHitTo;
        this.hspQueryFrame = hspQueryFrame;
        this.hspHitFrame = hspHitFrame;
        this.hspIdentity = hspIdentity;
        this.hspPositive = hspPositive;
        this.hspGaps = hspGaps;
        this.hspAlignLen = hspAlignLen;
        this.hspQseq = hspQseq;
        this.hspHseq = hspHseq;
        this.hspIdentityString = hspIdentityString;
    }
      private int hspQueryFrame;
      private int hspHitFrame;
      private int hspIdentity;
      private int hspPositive;
      private int hspGaps;
      private int hspAlignLen;
      private String hspQseq;
      private String hspHseq;
      private String hspIdentityString;
}
