package org.biojava.bio.symbol;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.io.CharacterTokenization;
import org.biojava.bio.seq.io.SymbolListCharSequence;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;


/**
 * <p>
 * A suffix tree is an efficient method for encoding the frequencies
 * of motifs in a sequence.  They are sometimes used to quickly screen
 * for similar sequences.  For instance, all motifs of length up to
 * 2 in the sequence <code>AAGT</code> could be encoded as:
 * </p>
 *
 * <pre>
 * root(4)
 * |
 * A(2)--------G(1)-----T(1)
 * |           |
 * A(1)--G(1)  T(1)
 * </pre>
 *<p>
 * It supports addition of elements both as String and SymbolList. They should
 * not be mixed together. The strings are also terminated internally, so it is
 * possible to go and add more than one string to the suffix tree.
 *</p>
 *<p>
 * Some more work need to be done on how data should be generated from this
 * class. If you need something that's not in there, please e-mail the list at
 * biojava-dev@biojava.org and I'll add it in there.
 *<\p>
 * @author Francois Pepin
 * @version 1.3
 */
public class UkkonenSuffixTree{

  public static final char DEFAULT_TERM_CHAR ='$';
  private char terminationChar;

  SuffixNode root;

  public static final int TO_A_LEAF = -1;
  private int e;

  private CharSequence sequences;


  /** Describes the rule that needs to be applied after walking down a tree. Put as a class variable because it can only return a single object (and I don't want to extend Node any further.
   * rule 1: ended up at a leaf.
   * rule 2: need to extend an internalNode.
   * rule 3: would split an edge.
   * rule 4: ended up in the middle of an edge.
   * rule 5: ended up at an InternalNode
   *
   * Production 5 counts as rule 4 when adding a sequence, but the rules are also used to when searching the tree.
   */
  private int rule;


  /** Initializes a new <code>UkkonenSuffixTree</code> instance.
   */
  public UkkonenSuffixTree(){
    terminationChar = DEFAULT_TERM_CHAR;
    root = new SimpleNode();
    e=0;
    sequences = "";
  }

  public UkkonenSuffixTree(String seqs){
    this();
    addSequence(seqs, "unnamed", false);
  }

  public UkkonenSuffixTree(FiniteAlphabet alpha){
    this();
  }

  public void addSymbolList(SymbolList list, String name, boolean doNotTerminate) throws IllegalSymbolException{
    if (!doNotTerminate)
      list = new TerminatedSymbolList(list);

    CharSequence seq =new SymbolListCharSequence(list);
    System.out.println("Adding symbol list " + seq.toString());
    //if (!doNotTerminate)

    //if(!doNotTerminate)
    //  seq = seq+terminationChar;
    addPreppedSequence(seq);
  }



  /**
   * Add a sequence into the tree. If there are more sequences, they should be separated by a terminationChar ($ by default). If none exist, it is assumed that there is only 1 continuous sequence to be added.
   * @param seq the sequence/sequences to be added into the tree.
   * @param doNotTerminate whether we should terminate the sequence if it's non-terminated.
   */
  public void addSequence(String seq, String name, boolean doNotTerminate){
    int i;
    int start, end;
    ArrayList toBeAdded = new ArrayList();
    Iterator iterator;
    String subseq;

    if (seq==null||seq.length()==0)
      return;

    //terminate the String if it's not terminated.
    if(!doNotTerminate&&seq.charAt(seq.length()-1)!=terminationChar)
      seq = seq+terminationChar;

    //count how many termination Chars in in.
    start =0;
    for (i=0;seq.indexOf(terminationChar, i)!=-1;i=seq.indexOf(terminationChar, i)+1){
      start=i;
      end=seq.indexOf(terminationChar, i);
      toBeAdded.add(seq.substring(start,end+1));
    }

    iterator = toBeAdded.iterator();
    i=0;
    while (iterator.hasNext()){
      subseq=(String)iterator.next();
      addPreppedSequence(subseq);
      i++;
    }
  }

  /** Add a single sequence into the tree.
   *
   * @param seq a <code>String</code> value
   */
  private void addPreppedSequence(CharSequence seq){
    int i, gammaStart;
    int j=0;
    SuffixNode oldNode=null, newNode;
    SuffixNode currentNode;
    boolean canLinkJump = false;

    //Puts i at the end of the previous sequences
    i = sequences.length();
    j=i;

    sequences= sequences.toString()+seq.toString();

    currentNode = root;
    //phase i
    for (; i<sequences.length();i++){
      //System.out.println("Phase "+i);

      e+=1;
      //extension j;
      for (;j<=i;j++){
	//System.out.println("extension "+j);

	//reset a couple of things...
	newNode = null;

	//find first node v at or above s[j-1,i] that is root or has a suffixLink
	while (currentNode!=root&&currentNode.suffixLink==null&&canLinkJump)
	  currentNode = currentNode.parent;

	if (root==currentNode){
          currentNode=jumpTo(root,sequences,j,i+1);
        }else{
          if (canLinkJump)
            currentNode = currentNode.suffixLink;
          gammaStart = j+getPathLength(currentNode);

          currentNode = jumpTo(currentNode,sequences,gammaStart,i+1);
	}

	if (rule==1)
	  addPositionToLeaf(j, currentNode);
	if (rule==2)
          doRule2(currentNode,i,j);
	if (rule==3){
	  newNode=doRule3(currentNode,i,j);
	  currentNode=newNode;
	}

	if (rule==1||rule==4||rule==5)
	  currentNode = currentNode.parent;

	if (oldNode!=null){
	  if (currentNode.isTerminal())
	    currentNode=currentNode.parent;

	  oldNode.suffixLink=currentNode;

	}
	oldNode=newNode;
	newNode=null;

        if (rule==1||rule==4||rule==5){
	  oldNode=null;
          canLinkJump=false;
          break;
	}else
          canLinkJump=true;


      }//for phase i
    }//for extension j
    finishAddition();
  }

  /** This method is used to walk down the tree, from a given node. The
   * <code>rule</code> variable can be used to check where the walk
   *  stopped. Note that rule 3 means that the string used to walk down the
   *  tree does not match (which is a bit different from the construction
   *  where rule 3 implies that only the last character doesn't match.
   *<p>
   *  The String is encoded as a substring of a given source. This is done to
   *  avoid replicating the string. To send walk down the string
   *  <code>x</code> from the root, one would call walkTo(root,x,0,x.length()).
   *
   * @param starting the root of the subtree we're walking down form.
   * @param source a superstring that contains the string we're using to
   * walking down. source.subtring(from,to) should give the string we're
   * walking down from.
   * @param from the start position (inclusive) of the target string in the
   * source.
   * @param to the end position (exclusive) of the target string in the node.
   * @return a <code>SuffixNode</code> that the walk stopped at. If the walk
   * stopped inside an edge. (check the rule variable to see where it stopped).
   */
  public SuffixNode walkTo(SuffixNode starting, String source, int from, int to){
    SuffixNode currentNode;
    SuffixNode arrivedAt;
    CharSequence edgeLabel;


    currentNode=starting;
    arrivedAt=starting;
    while (from<to){
      arrivedAt=(SuffixNode)currentNode.children.get(
                                                     new Character(source.charAt(from)));
      if (arrivedAt==null){
	from=to;
	arrivedAt=currentNode;
	rule=2;
	break;
      }

      edgeLabel = getEdgeLabel(arrivedAt);
      if (edgeLabel.length()>=to-from){
        if (edgeLabel.equals(source.substring(from,to))){
          //rule 1 or 5,
          if (arrivedAt.isTerminal())
            rule=1;
          else
            rule=5;
        }
        if (edgeLabel.subSequence(0, to-from).equals(source.substring(from,to)))
          rule=4;
        else
          rule=3;
        from=to;
      } else if (source.subSequence(from,from+edgeLabel.length())
                 .equals(edgeLabel)) {
        from+=edgeLabel.length();
        currentNode=arrivedAt;
      }

      else{
        rule=3;
        from=to;
      }
    }

    return arrivedAt;

  }




  /**
   * Just like walkTo, but faster when used during tree construction, as it
   * assumes that a mismatch can only occurs with the last character of the
   * target string.
   *
   * @param starting the root of the subtree we're walking down form.
   * @param source a superstring that contains the string we're using to
   * walking down. source.subtring(from,to) should give the string we're
   * walking down from.
   * @param from the start position (inclusive) of the target string in the
   * source.
   * @param to the end position (exclusive) of the target string in the node.
   * @return a <code>SuffixNode</code> that the walk stopped at. If the walk
   * stopped inside an edge. (check the rule variable to see where it
   * stopped).
   */
  public SuffixNode jumpTo(SuffixNode starting, CharSequence source, int from, int to){
    SuffixNode currentNode;
    SuffixNode arrivedAt;
    boolean canGoDown = true;
    int edgeLength;
    int original=from;
    SuffixNode originalNode=starting;


    currentNode=starting;
    arrivedAt=starting;

    rule=0;

    if (from==to){
      rule=5;
      return starting;
    }


    while (canGoDown){
      //    if (source.substring(from, to).equals("CAGCG"))
      //  System.out.println(to+" here to "+source.substring(from, to)+" "+(i++));

      if (currentNode.isTerminal()){
      	System.out.println("ARRGH! at "+source.subSequence(original, to)+
                           "("+from+","+original+","+to+
                           ") from "+getLabel(originalNode));
        //Something truly awful happened if this line is ever reached.
        //This bug should be dead, but it it came back from the dead a couple
        //of times already.
      }

      arrivedAt=(SuffixNode)currentNode.children.get(
                                                     new Character(source.charAt(from)));
      if (arrivedAt==null){
	canGoDown=false;
	arrivedAt=currentNode;
	rule=2;
	break;
      }

      edgeLength = getEdgeLength(arrivedAt);
      if (edgeLength>=to-from){
	int after = getPathEnd(arrivedAt)-getEdgeLength(arrivedAt)+to-from-1;
	if (sequences.charAt(after)==
	    source.charAt(to-1)){
	  if (getEdgeLength(arrivedAt)==to-from){
	    if (arrivedAt.isTerminal())
	      rule=1;
	    else
	      rule=5;
	  }
	  else
	    rule=4;
	}
	else
	  rule=3;
	canGoDown=false;
	break;
      }
      from+=edgeLength;
      currentNode=arrivedAt;

    }//while canGoDOwn

    return arrivedAt;
  }

  /******************************************************************
   * Tree navigation methods
   ******************************************************************/

  protected int getEdgeLength(SuffixNode child){
    int parentLength, childLength;
    SuffixNode parent;
    if (child==root)
      return 0;
    parent=child.parent;
    parentLength = getPathLength(parent);
    childLength = getPathLength(child);
    if (childLength-parentLength<=0){

      System.out.println("negative length "+(childLength-parentLength));

      System.out.println(getLabel(child)+","+getLabel(parent));
    }

    return childLength-parentLength;
  }

  protected CharSequence getEdgeLabel(SuffixNode child){
    return sequences.subSequence(
                               child.labelStart+
                               (getPathLength(child)-getEdgeLength(child)),
                               (child.labelEnd==TO_A_LEAF)?
                               e:child.labelEnd);
  }


  protected int getPathLength(SuffixNode node){
    return getPathEnd(node)-node.labelStart;
  }

  protected int getPathEnd(SuffixNode node){
    return node.labelEnd==TO_A_LEAF?e:node.labelEnd;
  }

  protected CharSequence getLabel(SuffixNode node){
    if (node==root)
      return "root";
    else
      return sequences.subSequence(
                                 node.labelStart,
                                 (node.labelEnd==TO_A_LEAF)?e:node.labelEnd).toString();
  }


  protected ArrayList getAllNodes(SuffixNode root, ArrayList list, boolean leavesOnly){
    Iterator iterator;
    if (list==null)
      list= new ArrayList();
    if (!leavesOnly||(leavesOnly&&root.isTerminal()))
      list.add(root);
    if (!root.isTerminal()){
      iterator = root.children.values().iterator();
      while (iterator.hasNext())
	list=getAllNodes((SuffixNode)iterator.next(), list, leavesOnly);
    }

    return list;
  }

  public void printTree(){
    ArrayList allNodes = getAllNodes(root, null, false);
    for (int i=0;i<allNodes.size();i++){
      SuffixNode node = (SuffixNode)allNodes.get(i);
      if (node==root)
        System.out.println("root");
      else
        System.out.println("node "+i+" label "+getLabel(node)+" attached to "+getLabel(node.parent));
    }
  }


  public SuffixNode getRoot(){return root;}

  /******************************************************************
   * End Tree Navigation Methods
   ******************************************************************/


  /******************************************************************
   * Tree modification methods
   ******************************************************************/
  private void addPositionToLeaf(int pos, SuffixNode leaf){
    int[] moreLabels;
    if (leaf.additionalLabels==null)
      leaf.additionalLabels = new int[]{pos};
    else{
      moreLabels = new int[leaf.additionalLabels.length+1];
      System.arraycopy(leaf.additionalLabels, 0, moreLabels, 0, leaf.additionalLabels.length);
      moreLabels[moreLabels.length-1]=pos;
      leaf.additionalLabels=moreLabels;
    }

  }

  private void doRule2(SuffixNode parent, int splittingPos, int suffixStart){
    //int number = getAllNodes(root, null, false).size();
    SuffixNode leaf = new SimpleNode (parent, suffixStart);

    parent.children.put(new Character(sequences.charAt(splittingPos)), leaf);
    //System.out.println("rule 2: "+sequences.charAt(splittingPos)+" from "+getLabel(parent)+ " to "+getLabel(leaf));

  }

  private SuffixNode doRule3(SuffixNode child, int splittingPos, int suffixStart){
    //      return toBeSplit.splitEdge(endOfSubSeq, sequences.charAt(endOfSubSeq),
    //			       toBeSplit.getStart()+endOfSubSeq-rule3Position,
    //			       suffixStart);
    //int number = getAllNodes(root, null, false).size();
    SuffixNode parent = child.parent;
    SuffixNode middle= new SimpleNode(parent,suffixStart,splittingPos);
    Character x=new Character(
                              sequences.charAt(child.labelStart+getPathLength(child)-getEdgeLength(child)));

    //System.out.println(parent.children.get(x)==child);

    Character y=new Character(sequences.charAt(
                                               child.labelStart
                                               +getPathLength(child)-getEdgeLength(child)
                                               +getEdgeLength(middle)
                                               ));

    parent.children.remove(x);
    parent.children.put(x,middle);

    middle.children.put(y,child);
    child.parent=middle;
    //System.out.println("rule 3: "+sequences.charAt(splittingPos)+" between "+getLabel(parent)+" and "+getLabel(child) + " Addition made:"+(number==getAllNodes(root, null,false).size()-1));
    doRule2(middle,splittingPos,suffixStart);
    return middle;
  }

  private void finishAddition(){
    SuffixNode leaf;
    ArrayList leaves = getAllNodes(root, null, true);
    for (int i=0;i<leaves.size();i++){
      leaf = (SuffixNode)leaves.get(i);
      if (leaf.labelEnd==TO_A_LEAF)
	leaf.labelEnd=e;
    }

  }
  /******************************************************************
   * end Tree modification methods
   ******************************************************************/

  public static abstract class SuffixNode {

    static final int A_LEAF=-1;
    SuffixNode parent;
    SuffixNode suffixLink;
    int labelStart, labelEnd;
    HashMap children;
    int[] additionalLabels;

    /**
       * Determine is this node is terminal (has no children).
       *<p>
       * Note that this only happens at the terminated node (if the sequences
       * have been terminated.
       *
       * @return <code>true</code> if and only if this node has no children.
       */

    abstract public boolean isTerminal();

    /**
     * Determine if this node has a child corresponding to a given character
     *
     * @param i the first <code>Character</code> of the edge coming down this node.
     * @return <code>true</code> if the node has a child going down from that character,
     * false otherwise
     */
    abstract public boolean hasChild(Character i);

    /** Gets the child of of a node that is coming down from that particular
     * node. It returns null if no child exists or if no child exists starting
     * on that particular character.
     *
     * @param i the first <code>Character</code> of the edge coming down this node
     * @return the appropriate child <code>SuffixNode</code>, or null if it
     * doesn't exists.
     */
    abstract SuffixNode getChild(Character i);
    //abstract void addChild(SuffixTree tree, int i, SuffixNode n);

    /**
     * Returns the parent of this node, null if it's the root.
     *
     * @return the parent of this node, null if it's the root.
     */
    abstract SuffixNode getParent();

  }


  class SimpleNode extends SuffixNode{

    //static final int A_LEAF=-1;
    //SuffixNode parent;
    //SuffixNode suffixLink;
    //int labelStart, labelEnd;
    //HashMap children;
    //int[] additionalLabels;

    /** Creates a root
     */
    public SimpleNode(){
      parent=null;
      suffixLink=null;
      labelStart=0;
      labelEnd=0;
      children=new HashMap();
      additionalLabels=null;
    }

    /** creates a leaf
     * @param parent the parent node
     * @param position the starting value of the suffix
     */
    public SimpleNode(SuffixNode parent, int position){
      this();
      this.parent=parent;
      labelStart=position;
      labelEnd = A_LEAF;
      children=null;
      //checkParent(this);
    }

    /** creates an internal node
     * @param parent the parent of this node
     * @param labelStart the starting point of the path label
     * @param labelStop the ending point of the path label
     */
    public SimpleNode(SuffixNode parent, int labelStart, int labelStop){
      this();
      this.parent=parent;
      this.labelStart=labelStart;
      this.labelEnd=labelStop;
      //checkParent(this);
    }


    public boolean isTerminal(){return children==null;}
    public boolean hasChild(Character x){return getChild(x)!=null;}
    public SuffixNode getChild(Character x){
      return (children==null)?null:(SuffixNode)children.get(x);
    }
    public SuffixNode getParent(){return parent;}
  }

  public boolean subStringExists(String str)
  {
    walkTo(root, str, 0, str.length());
    return (rule==1||rule==4||rule==5);
  }

  /**
   * Stupid little wrapper to put a termination symbol at the end. A sublist
   * that doesn't include the termination symbol will have the same alphabet as
   * the original, one that does will have an alphabet that includes the
   * termination symbol.
   *
   * Now includes some ugly hacks by Thomas Down to make ambiguities work
   * nicer.
   *
   */
  private class TerminatedSymbolList implements SymbolList
  {
    private SymbolList unterminated;
    final Symbol TERMINATION_SYMBOL;
    private AbstractAlphabet alpha;
    private Map translationTable = new HashMap();

    public TerminatedSymbolList(SymbolList unterminated)
    {
      this.unterminated=unterminated;
      TERMINATION_SYMBOL = AlphabetManager.createSymbol("Termination");
      FiniteAlphabet oldAlphabet = (FiniteAlphabet)unterminated.getAlphabet();
      // Set set = AlphabetManager.getAllSymbols(oldAlphabet);
      Set set = new HashSet();
      for (Iterator i = oldAlphabet.iterator(); i.hasNext(); ) {
          set.add(i.next());
      }
      set.add(TERMINATION_SYMBOL);
      alpha = new SimpleAlphabet(set);
      CharacterTokenization tokenizer =new CharacterTokenization(alpha, true);
      tokenizer.bindSymbol(TERMINATION_SYMBOL, DEFAULT_TERM_CHAR);
      SymbolTokenization sToke;
      try{
        sToke = oldAlphabet.getTokenization("token");
      }
      catch (BioException be){
          throw new BioError("Internal error: failed to get SymbolTokenization for SymbolList alphabet", be);
      }
      if (sToke.getTokenType() != SymbolTokenization.CHARACTER)
        throw new IllegalArgumentException("Only FiniteAlphabets using a char token are supported by UkkonenSuffixTree");
      try{
        for (Iterator i= AlphabetManager.getAllSymbols(oldAlphabet).iterator();i.hasNext();){
          Symbol oldSymbol = (Symbol) i.next();
          Symbol newSymbol;
          if (oldSymbol instanceof AtomicSymbol) {
              // The atomic symbols are identical between the two alphabets
              newSymbol = oldSymbol;
          } else {
              // Port ambiguous symbols across
              Set s = new HashSet();
              for (Iterator si = ((FiniteAlphabet) oldSymbol.getMatches()).iterator(); si.hasNext(); ) {
                  s.add(si.next());
              }
              newSymbol = alpha.getAmbiguity(s);
          }
          //takes first char of String, should work because we checked
          //getTokenType above

          char c = sToke.tokenizeSymbol(oldSymbol).charAt(0);
          // System.err.println("Binding " + c);
          tokenizer.bindSymbol(newSymbol, c);
          translationTable.put(oldSymbol, newSymbol);
        }
        //This is really hacky, ambiguous symbols containing
        //TERMINATION_SYMBOL are impossible at this point, so I just define
        //the Tokenization to treat them as TERMINATION_SYMBOL so that the
        //code that likes to loop through doesn't go titty up.
        for (Iterator i= AlphabetManager.getAllSymbols(alpha).iterator();i.hasNext();){
          Symbol s = (Symbol) i.next();
          Alphabet mathes = s.getMatches();
          if (mathes.contains(TERMINATION_SYMBOL))
              tokenizer.bindSymbol(s,DEFAULT_TERM_CHAR);
        }

      }catch(IllegalSymbolException ise){
        throw new AssertionFailure("Assertion Failure: This alphabet has been custom made so this doesn't happen",ise);
      }

      alpha.putTokenization("token", tokenizer);
    }


      // Implementation of org.biojava.utils.Changeable

    public void addChangeListener(ChangeListener changeListener, ChangeType changeType) {
      unterminated.addChangeListener(changeListener,changeType);
    }

    public void addChangeListener(ChangeListener changeListener) {
      unterminated.addChangeListener(changeListener);
    }

    public void removeChangeListener(ChangeListener changeListener, ChangeType changeType) {
      unterminated.removeChangeListener(changeListener,changeType);
    }

    public void removeChangeListener(ChangeListener changeListener) {
      unterminated.removeChangeListener(changeListener);
    }

    public boolean isUnchanging(ChangeType changeType) {
      return unterminated.isUnchanging(changeType);
    }

    // Implementation of org.biojava.bio.symbol.SymbolList

    public int length() {
      return unterminated.length()+1;
    }

    public Iterator iterator() {
      return unterminated.iterator();
    }

    public SymbolList subList(int n, int n1) throws IndexOutOfBoundsException {
      List list;
      if (n1!=unterminated.length()+1)
        return unterminated.subList(n,n1);
      else{
        list = unterminated.subList(n,n1-1).toList();
        list.add(TERMINATION_SYMBOL);
        try{
          return new SimpleSymbolList(getAlphabet(),list);
        }
        catch(IllegalSymbolException e){
          throw new AssertionFailure("Assertion Failure: This alphabet was created just so it doesn't do this",e);
        }
      }
    }

    public Alphabet getAlphabet() {
      return alpha;
    }

    public Symbol symbolAt(int n) throws IndexOutOfBoundsException {
      if (n!=length())
        return (Symbol) translationTable.get(unterminated.symbolAt(n));
      else
        return TERMINATION_SYMBOL;
    }

    public List toList() {
      List answer = unterminated.toList();
      answer.add(TERMINATION_SYMBOL);
      return answer;
    }

    public String seqString() {
      try{
        SymbolTokenization toke = getAlphabet().getTokenization("token");
        return unterminated.seqString()+toke.tokenizeSymbol(TERMINATION_SYMBOL);
      } catch (BioException ex) {
        throw new BioRuntimeException("Couldn't tokenize sequence", ex);
      }
    }

      public String subStr(int n, int n1) throws IndexOutOfBoundsException {
      return subList(n, n1).seqString();
    }

    /**
     * Describe <code>edit</code> method here.
     *
     * @param edit an <code>Edit</code> value
     * @exception IndexOutOfBoundsException if an error occurs
     * @exception IllegalAlphabetException if an error occurs
     * @exception ChangeVetoException if an error occurs
     */
    public void edit(Edit edit) throws IndexOutOfBoundsException, IllegalAlphabetException, ChangeVetoException {
      throw new ChangeVetoException("TerminatedSymbolList is immutable");
    }
  }
}
