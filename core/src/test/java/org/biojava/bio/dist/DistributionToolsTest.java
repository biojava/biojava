package org.biojava.bio.dist;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import junit.framework.TestCase;

import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Symbol;

/**
 * Tests that methods from DistributionTools work as advertised.
 *
 * @author Mark Schreiber
 * @since 1.3
 */
public class DistributionToolsTest extends TestCase {

  private Alignment a;
  private Distribution random;
  private String[] sa;
  private Map map;

  public DistributionToolsTest(String name) {
    super(name);
  }

  protected void setUp() {
    try{
      sa = new String[]{"CA-TGGG","AATTGGC","AATTGGG",
                        "AATTGGC","AA-TGGG","AATTGGC",
                        "AATTGGG","AATTGGC","AATTGGG",
                        "AATTGGC"};

      map = new HashMap(sa.length);
      for (int i = 0; i < sa.length; i++) {
         map.put(new Integer(i), DNATools.createDNA(sa[i]));
      }

      a = new SimpleAlignment(map);

      random = new SimpleDistribution(DNATools.getDNA());
      DistributionTools.randomizeDistribution(random);

    }catch(Exception e){
      e.printStackTrace();
    }
  }

  protected void tearDown(){
    sa = null;
    map = null;
    a = null;
    random = null;
  }


  public void testDistOverAlignment() {
    try{
      Distribution[] d = DistributionTools.distOverAlignment(a,false);
      Distribution[] d2 = DistributionTools.distOverAlignment(a,false,10.0);
      Distribution[] d3 = DistributionTools.distOverAlignment(a,true);

      assertTrue(d[0].getWeight(DNATools.a()) == 0.9);
      assertTrue(d[0].getWeight(DNATools.c()) == 0.1);
      assertTrue(d[0].getWeight(DNATools.g()) == 0.0);
      assertTrue(d[0].getWeight(DNATools.t()) == 0.0);
      assertTrue(d[0].getWeight(AlphabetManager.getGapSymbol()) == 0.0);

      assertTrue(d2[1].getWeight(DNATools.a()) == 0.625);
      assertTrue(d2[1].getWeight(DNATools.c()) == 0.125);
      assertTrue(d2[1].getWeight(DNATools.g()) == 0.125);
      assertTrue(d2[1].getWeight(DNATools.t()) == 0.125);
      assertTrue(d2[1].getWeight(AlphabetManager.getGapSymbol()) == 0.0);

      assertEquals( 0.0 , d3[2].getWeight(DNATools.a()) ,0.000001);
      assertEquals( 0.0 , d3[2].getWeight(DNATools.c()) ,0.000001);
      assertEquals( 0.0 , d3[2].getWeight(DNATools.g()) ,0.000001);
      assertEquals( 0.8 , d3[2].getWeight(DNATools.t()), 0.000001);
      assertEquals( 0.2 , d3[2].getWeight(AlphabetManager.getGapSymbol()),
                    0.000001);


    }catch(Exception e){
      e.printStackTrace();
    }
  }

  /**
   * This method tests DistributionTools.shannonEntropy(), DistributionTools.totalEntropy()
   * and DistributionTools.bitsOfInformation().
   */
  public void testInformationContent(){
    Distribution d = new UniformDistribution(DNATools.getDNA());
    assertTrue(DistributionTools.bitsOfInformation(d) == 0.0);
    assertTrue(DistributionTools.totalEntropy(d) == 2.0);

    Map map = DistributionTools.shannonEntropy(d,2.0);
    for(Iterator i = DNATools.getDNA().iterator(); i.hasNext();){
      Symbol s = (Symbol)i.next();
      double ent = ((Double)map.get(s)).doubleValue();
      assertTrue(ent == 2.0);
    }
  }

  public void testGenerateSequence(){
    Sequence seq = DistributionTools.generateSequence("seq",random,1000);
    assertTrue(seq.length() == 1000);
    //this will kick in the PackedSymbolListFactory
    Sequence seq2 = DistributionTools.generateSequence("seq", random, 30000);
    assertTrue(seq2.length() == 30000);
  }

}
