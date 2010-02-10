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

package org.biojava.bio.gui;

import java.util.ArrayList;
import java.util.Iterator;

import javax.swing.JTree;
import javax.swing.event.TreeModelListener;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreePath;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.SequenceDB;

/**
 * FeatureTree is GUI tree to display the features and annotations
 * of the sequences in a <code>SequenceDB</code> Nested Features are
 * displayed as expandable leaves.
 *
 * <p>Copyright:    Copyright (c) 2002</p>
 * <p>Company:      AgResearch</p>
 *
 * @author Mark Schreiber
 * @version 1.0
 */

public class FeatureTree extends JTree{
  private String root = "DB";
  private ArrayList seqs = new ArrayList();

  /**
   * Create a new FeatureTree
   */
  public FeatureTree() {
    super();
    this.setModel(new FeatureModel(root));
  }

  /**
   * Use this method to provide the sequences for the tree to work with.
   * @param db A database of <CODE>Sequence</CODE>s to display
   * @throws org.biojava.bio.BioException if information cannot be retrieved from <CODE>db</CODE>
   */
  public void setSequenceDB(SequenceDB db) throws BioException{
    for(SequenceIterator i = db.sequenceIterator(); i.hasNext();){
      seqs.add(i.nextSequence());
    }
    this.expandRow(0);
    this.repaint();
  }

  /**
   * Labels <code>Sequence</code> objects with their name, <code>Annotations</code> with
   * the tag Annotations, <code>Features</code> with the tag Features and other objects
   * with the <code>toString</code> value.
   */
  public String convertValueToText(Object value, boolean selected,
                                     boolean expanded, boolean leaf, int row,
                                     boolean hasFocus) {
      if(value instanceof Sequence) return ((Sequence)value).getName();
      if(value instanceof Annotation) return "Annotations";
      if(value instanceof Feature) return value.toString();
      if(value instanceof FeatureHolder) return "Features";
      else if(value != null) return value.toString();
      return "";
    }

  /**
   * The model used by the FeatureTree
   */
  class FeatureModel implements TreeModel {
    protected String root;

    public FeatureModel(String root){
      this.root = root;
    }

    public Object getRoot(){ return root;}

    public boolean isLeaf(Object node){
      if(node.equals(getRoot())){
        if(seqs.size()==0) return true;
        else return false;
      }
      if(node instanceof FeatureHolder){
        return false;
      }
      if(node.equals(Annotation.EMPTY_ANNOTATION)){
        return true;
      }
      if(node instanceof String) return true;
      return false;
    }

    public int getChildCount(Object parent){
      if(parent.equals(this.getRoot())) return seqs.size();
      else if(parent instanceof Sequence){
        return 3; //annotation + feature holder + SymbolList
      }
      if(parent instanceof Annotation){
        return ((Annotation)parent).keys().size();
      }
      if(parent instanceof Feature){
        return 3;//for annotation + Sequence + Nested Features
      }
      if(parent instanceof FeatureHolder){
        return ((FeatureHolder)parent).countFeatures();
      }

      return 0;
    }

    public Object getChild(Object parent, int index){
     if(parent.equals(getRoot())){
       return ((Sequence)seqs.get(index));
     }else if(parent instanceof Sequence){
       if(index == 0) return ((Sequence)parent).getAnnotation();
       if(index == 1) return ((Sequence)parent).filter(FeatureFilter.all,false);
       else return ((Sequence)parent).seqString();
     }else if(parent instanceof Feature){
       if(index == 0)return ((Feature)parent).getAnnotation();
       if(index == 1)return ((Feature)parent).filter(FeatureFilter.all,false);
       else return ((Feature)parent).getSymbols().seqString();

     }else if(parent instanceof FeatureHolder){
       ArrayList al = new ArrayList();
       for(Iterator i = ((FeatureHolder)parent).features(); i.hasNext();){
         al.add(i.next());
       }

       Feature f = (Feature)al.get(index);
       return f;
     }else{//parent must be Annotation
       ArrayList al = new ArrayList();
       for(Iterator i = ((Annotation)parent).keys().iterator(); i.hasNext();){
         Object key = i.next();
         Object value = ((Annotation)parent).getProperty(key);
         al.add(key.toString()+" : "+value.toString());
       }
       return al.get(index);
     }
    }

    public int getIndexOfChild(Object parent, Object child){
     if(parent.equals(getRoot())){
       for(int i = 0; i < seqs.size(); i++){
         if(seqs.get(i).equals(child)) return i;
       }
     }else if(parent instanceof Sequence || parent instanceof Feature){
       if(child instanceof Annotation)return 0;
       if(child instanceof FeatureHolder)return 1;
       else return 2;
     }else if(parent instanceof FeatureHolder){
       ArrayList al = new ArrayList();
       for(Iterator i = ((FeatureHolder)parent).features(); i.hasNext();){
         al.add(i.next());
       }
       for(int i = 0; i < al.size(); i++){
         if(al.get(i).equals(child)) return i;
       }
     }else if(parent instanceof Annotation){
              ArrayList al = new ArrayList();
       for(Iterator i = ((Annotation)parent).keys().iterator(); i.hasNext();){
         Object key = i.next();
         Object value = ((Annotation)parent).getProperty(key);
         al.add(key.toString()+" : "+value.toString());
       }
       for (int i = 0; i < al.size(); i++) {
         if(child.equals(al.get(i))) return i;
       }

     }
      return -1;
    }

    public void valueForPathChanged(TreePath path, Object newValue){}
    public void removeTreeModelListener(TreeModelListener l){}
    public void addTreeModelListener(TreeModelListener l){}
  }
}
