/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 */


package org.biojava.bio;


import junit.framework.TestCase;
import junit.framework.TestSuite;

/** List of tests for for MergeAnnotationTree.
 * @author Francois Pepin
 * @version 1.4
 */
public class MergeAnnotationTest extends TestCase {

  Annotation anno1;
  Annotation anno2;

  MergeAnnotation merged;

  public static void main(String args[]){
    junit.textui.TestRunner.run(new TestSuite(MergeAnnotationTest.class));
  }
  
  protected void setUp() throws Exception
  {
    anno1 = new SimpleAnnotation();
    anno1.setProperty("name", "michelle");
    anno1.setProperty("age", "27");
        
    anno2 = new SmallAnnotation();
    anno2.setProperty("work", "mcb");
    anno2.setProperty("favorite food", "?");
    anno2.setProperty("cat", "bretelle&trouser");
    
    merged = new MergeAnnotation();
    merged.addAnnotation(anno1);
    merged.addAnnotation(anno2);
  }

  public void testAnnotation() {
    assertEquals(2,merged.getAnnotations().size());
  }
  
  
  public void testKeys()
  {
    assertEquals(5, merged.keys().size());
  }
  
}

      
