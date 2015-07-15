/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.core.search.io;

import java.util.ArrayList;
import java.util.Iterator;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 *
 * @author pavanpa
 */
public  class Hit implements Iterable<Hsp>{
    private final int hitNum;
    private final String hitId;
    private final String hitDef;
    private final String hitAccession;
    private final int hitLen;
    private final ArrayList<Hsp> hsps;
    private Sequence hitSequence;
    
    

    public Hit(int hitNum, String hitId, String hitDef, String hitAccession, int hitLen, ArrayList<Hsp> hsps, Sequence hitSequence) {
        this.hitNum = hitNum;
        this.hitId = hitId;
        this.hitDef = hitDef;
        this.hitAccession = hitAccession;
        this.hitLen = hitLen;
        this.hsps = hsps;
        this.hitSequence = hitSequence;
    }

    public int getHitNum() {
        return hitNum;
    }

    public String getHitId() {
        return hitId;
    }

    public String getHitDef() {
        return hitDef;
    }

    public String getHitAccession() {
        return hitAccession;
    }

    public int getHitLen() {
        return hitLen;
    }
    
    /**
     * returns the reference to the original and whole sequence hit in the database.
     * Available only if the ResultFactory implements setHitReferences and
     * it was used before the parsing with SearchIO
     * @return Sequence object
     */
    public Sequence getHitSequence() {
        return hitSequence;
    }
    
    @Override
    public Iterator<Hsp> iterator() {
        return new Iterator<Hsp>() {
            int current = 0;
            @Override
            public boolean hasNext() {
                return current < hsps.size();
            }

            @Override
            public Hsp next() {
                return hsps.get(current++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException("The remove operation is not supported by this iterator");
            }
        };
    }
}
