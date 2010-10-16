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
 *
 */
package org.biojava.utils.bytecode;

import java.util.*;

/**
 * A list of Instructions and/or other CodeGenerator objects.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public class InstructionVector implements CodeGenerator {
  private final List instructions;
  private final Label startLabel;
  private final Label endLabel;
  
  private StatsCache statsCache; // gets set to null on edits
  
  {
    instructions = new ArrayList();
    startLabel = new Label();
    endLabel = new Label();
    statsCache = null;
  }
  
  public void add(CodeGenerator g) {
    statsCache = null;
    instructions.add(g);
  }
  
  public int size() {
    return instructions.size();
  }
  
  public void add(int pos, CodeGenerator g) {
    statsCache = null;
    instructions.add(pos, g);
  }
  
  public void remove(int pos) {
    statsCache = null;
    instructions.remove(pos);
  }
  
  public CodeGenerator generatorAt(int pos) {
    return (CodeGenerator) instructions.get(pos);
  }
  
  public Label getStartLabel() {
    return startLabel;
  }
  
  public Label getEndLabel() {
    return endLabel;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    CodeContext subctx = ctx.subContext();
    subctx.open();
    subctx.markLabel(startLabel);
    for (Iterator i = instructions.iterator(); i.hasNext(); ) {
      CodeGenerator cg = (CodeGenerator) i.next();
      cg.writeCode(subctx);
    }
    subctx.markLabel(endLabel);
    subctx.close(); // Wrap up the subcontt
  }
  
  public int stackDepth() {
    StatsCache statsCache = getStatsCache();
    
    return statsCache.depth;
  }
  
  public int stackDelta() {
    StatsCache statsCache = getStatsCache();
    
    return statsCache.delta;
  }
  
  private StatsCache getStatsCache() {
    if(statsCache == null) {
      int depth = 0;
      int delta = 0;
      
      for(Iterator i = instructions.iterator(); i.hasNext(); ) {
        CodeGenerator cg = (CodeGenerator) i.next();
        int dp = cg.stackDepth();
        int dl = cg.stackDelta();
        
        dp += delta;
        delta += dl;

        depth = Math.max(depth, dp);
      }
      
      statsCache = new StatsCache(depth, delta);
    }
    
    return statsCache;
  }
  
  private static class StatsCache {
    public final int depth;
    public final int delta;
    
    public StatsCache(int depth, int delta) {
      this.depth = depth;
      this.delta = delta;
    }
  }
}
