/*
 * DPSerializationTest.java
 * JUnit based test
 *
 * Created on September 27, 2007, 10:17 AM
 */

package org.biojava.bio.dp;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import junit.framework.*;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.onehead.SingleDP;
import org.biojava.bio.seq.DNATools;

/**
 * Test of fix for bug #2359
 * @author Mark Schreiber
 */
public class DPSerializationTest extends TestCase {
    
    public DPSerializationTest(String testName) {
        super(testName);
    }
    
    public void testDPSerialization() throws Exception{
        MarkovModel model = new ProfileHMM(DNATools.getDNA(),5,DistributionFactory.DEFAULT, DistributionFactory.DEFAULT, "test");
        DP singleDP = new SingleDP(model);
        
        MarkovModel model2 = new SimpleMarkovModel(2, DNATools.getDNA(), "pair");
        DP pairDP = DPFactory.DEFAULT.createDP(model2);
        
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        ObjectOutputStream oos = new ObjectOutputStream(baos);
        oos.writeObject(singleDP);
        oos.writeObject(pairDP);
        oos.flush();
        oos.close();
        
        
        ObjectInputStream ois = new ObjectInputStream(new ByteArrayInputStream(baos.toByteArray()));
        DP dp = (DP)ois.readObject();
        assertEquals(singleDP.getClass(), dp.getClass());
        
        DP dp2 = (DP)ois.readObject();
        assertEquals(pairDP.getClass(), dp2.getClass());
    }
}
