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
package org.biojavax.bio.phylo.io.nexus;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import java.util.Vector;

import org.jgrapht.WeightedGraph;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

import org.biojava.bio.seq.io.ParseException;


/**
 * Represents Nexus trees blocks.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @author Tiago Antao
 * @since 1.6
 */
public class TreesBlock extends NexusBlock.Abstract {

	/**
	 * A constant representing the name of Trees blocks.
	 */
	public static final String TREES_BLOCK = "TREES";

	private Map translations = new LinkedHashMap();

	private List comments = new ArrayList();

	private Map trees = new LinkedHashMap();

	private WeightedGraph<String, DefaultWeightedEdge> weighted =  new SimpleWeightedGraph<String, DefaultWeightedEdge>(DefaultWeightedEdge.class);

    private String topNode = null;

	private String nodePrefix = "p";

	private int pValue = 0;
	private Vector<String> uuids;

	/**
	 * A simple representation of a Newick tree as a single string.
	 */
	public static class NewickTreeString {
		private String rootType;

		private String treeString;

		private boolean starred;

		/**
		 * Make the tree (un)rooted.
		 * 
		 * @param rootType
		 *            'U' for unrooted, 'R' for rooted, <tt>null</tt> for
		 *            unsure.
		 */
		public void setRootType(final String rootType) {
			this.rootType = rootType;
		}

		/**
		 * Set the Newick string describing the tree.
		 */
		public void setTreeString(final String treeString) {
			this.treeString = treeString;
		}

		/**
		 * Sets whether this tree has a star before it's name.
		 * 
		 * @param starred
		 *            <tt>true</tt> if it has one.
		 */
		public void setStarred(boolean starred) {
			this.starred = starred;
		}

		/**
		 * Tests whether this tree has a star before it's name.
		 * 
		 * @return starred <tt>true</tt> if it has one.
		 */
		public boolean isStarred() {
			return this.starred;
		}

		/**
		 * See if the tree is rooted.
		 * 
		 * @return 'U' for unrooted, 'R' for rooted, <tt>null</tt> for unsure.
		 */
		public String getRootType() {
			return this.rootType;
		}

		/**
		 * Get the Newick string describing the tree.
		 * 
		 * @return the tree string.
		 */
		public String getTreeString() {
			return this.treeString;
		}
	}

	/**
	 * Delegates to NexusBlock.Abstract constructor using TreesBlock.TREES_BLOCK
	 * as the name.
	 */
	public TreesBlock() {
		super(TreesBlock.TREES_BLOCK);
	}

	/**
	 * Add a translation.
	 * 
	 * @param label
	 *            the label to add.
	 * @param taxa
	 *            the taxa name this label will represent.
	 */
	public void addTranslation(final String label, final String taxa) {
		this.translations.put(label, taxa);
	}

	/**
	 * Removes the given translation.
	 * 
	 * @param label
	 *            the label to remove.
	 */
	public void removeTranslation(final String label) {
		this.translations.remove(label);
	}

	/**
	 * Checks to see if we contain the given translation.
	 * 
	 * @param label
	 *            the label to check for.
	 * @return <tt>true</tt> if we already contain it.
	 */
	public boolean containsTranslation(final String label) {
		return this.translations.containsKey(label);
	}

	/**
	 * Get the translations added so far.
	 * 
	 * @return the translations added so far.
	 */
	public Map getTranslations() {
		return this.translations;
	}

	/**
	 * Adds a tree.
	 * 
	 * @param label
	 *            the label to give the tree.
	 * @param tree
	 *            the tree to add.
	 */
	public void addTree(final String label, final NewickTreeString tree) {
		this.trees.put(label, tree);
	}

	/**
	 * Removes a tree.
	 * 
	 * @param label
	 *            the label to remove.
	 */
	public void removeTree(final String label) {
		this.trees.remove(label);
	}

	/**
	 * Checks to see if we contain the given tree.
	 * 
	 * @param label
	 *            the label to check for.
	 * @return <tt>true</tt> if we already contain it.
	 */
	public boolean containsTree(final String label) {
		return this.trees.containsKey(label);
	}

	/**
	 * Returns all trees.
	 * 
	 * @return all the selected trees.
	 */
	public Map getTrees() {
		return this.trees;
	}

	/**
	 * Returns a tree for given label
	 * @param label
             * 	 the label to select.
	 *
	 * @return selected tree.
             */
     	public Object getTree(final String label) {
		return this.trees.get(label);
	}

	public String getTreeText(UndirectedGraph<String, DefaultEdge> treegraph, String node, String parent){
		Set<DefaultEdge> edges = treegraph.edgesOf(node);
		int numOffspring = 0;
		StringBuffer nodeText = new StringBuffer("");
		for(DefaultEdge e : edges) {
			String child;
			if (treegraph.getEdgeSource(e).equals(node)) {
				child = treegraph.getEdgeTarget(e);
			}
			else {
				child = treegraph.getEdgeSource(e);
			}
			if (child.equals(parent)) {
				break; //We dont want to go up the tree...
			}
			if (numOffspring>0) {
				nodeText.append(", ");
			}
			nodeText.append(getTreeText(treegraph, child, node));
			numOffspring++;
		}
		if (numOffspring==0) {
			return node;
		}
		else {
			return "(" + nodeText + ")";
		}
	}


	

	/**
	 * Add a tree, converting weighted graph (JGraphT) to NewickString
	 *
	 * This will assume a (arbitrary) root node using the old convention
	 * of labeling intermediate nodes as p*.
	 *
	 *
	 * @deprecated
	 * @param label
	 * 		  the label to add
	 *
	 * @param treegraph
	 * 		  the treegraph to convert.
     */      
	public void addTree(final String label, WeightedGraph<String, DefaultWeightedEdge> treegraph) {
		addTree(label, treegraph, "p0");
	}

	public String getTreeText(WeightedGraph<String, DefaultWeightedEdge> treegraph, String node, String parent){
		Set<DefaultWeightedEdge> edges = treegraph.edgesOf(node);
		int numOffspring = 0;
		StringBuffer nodeText = new StringBuffer("");
		for(DefaultWeightedEdge e : edges) {
			String child;
			if (treegraph.getEdgeSource(e).equals(node)) {
				child = treegraph.getEdgeTarget(e);
			}
			else {
				child = treegraph.getEdgeSource(e);
			}
			if (child.equals(parent)) {
				break; //We dont want to go up the tree...
			}
			if (numOffspring>0) {
				nodeText.append(", ");
			}
			nodeText.append(getTreeText(treegraph, child, node));
			nodeText.append(":"+treegraph.getEdgeWeight(e));
			numOffspring++;
		}
		if (numOffspring==0) {
			return node;
		}
		else {
			return "(" + nodeText + ")";
		}
	}

	/**
	 * Add a tree, converting weighted graph (JGraphT) to NewickString.
	 *
	 * @param label
	 * 		  the label to add
	 *
	 * @param treegraph
	 * 		  the treegraph to convert.
	 *
	 * @param topLabel
	 *        the label of the top (root if rooted tree) node.
     */
	public void addTree(final String label,
			WeightedGraph<String, DefaultWeightedEdge> treegraph, String topLabel) {
		final NewickTreeString tree = new NewickTreeString();
		String temp;

		for(String vertex : treegraph.vertexSet()) {
			if (vertex.equals(topLabel)) {
				topLabel = vertex; //String equality is not enough, has to be the same object
			}
		}
		temp = getTreeText(treegraph, topLabel, null);
		tree.setTreeString(temp);
		this.trees.put(label, tree);
	}	


	/**
	 * Renames a vertex of the weighted graph.
	 *
	 * @param oldName old name.
	 * @param newName new name.
	 */
	private void renameVertex(String oldName, String newName) {
		Set<DefaultWeightedEdge> oldEdges = this.weighted.edgesOf(oldName);

		this.weighted.addVertex(newName);
		for (DefaultWeightedEdge e : oldEdges) {
			DefaultWeightedEdge tempEdge;
			if (this.weighted.getEdgeSource(e).equals(oldName)) {
				tempEdge = this.weighted.addEdge(newName,
						this.weighted.getEdgeTarget(e));
			}
			else {
				tempEdge = this.weighted.addEdge(this.weighted.getEdgeSource(e),
						newName);
			}
			this.weighted.setEdgeWeight(tempEdge,this.weighted.getEdgeWeight(e));
			//this.unweighted.removeEdge(e); Not needed
		}
		this.weighted.removeVertex(oldName);
	}



	/**
	 * Tokenizes a string representing a newick tree.
	 * 
	 * Simple tokenizing, removing spaces and so on.
	 * 
	 * @param text the string representing the tree
	 * @return An array of tokens
	 */
	Vector<String> tokenize(String text) {
		Vector<String> toks = new Vector<String>();
		while (!text.equals("")) {
			text = text.trim();
			String delims = "(),: \t";
			String c = text.substring(0, 1);
			if (delims.contains(c)) {
				toks.add(c);
				text = text.substring(1);
			}
			else {
				String currTok = "";
				int pos = 0;
				while(!delims.contains(c)) {
					//StringBuffer here might be faster....
					currTok += c;
					if (text.length()>1) {
						text = text.substring(1);
						c = text.substring(0, 1);
					}
					else {
						toks.add(currTok);
						return toks;
					}
				}
				toks.add(currTok);
			}
		}

		return toks;
	}

	/*
	 * Checks if the graph has a name as vertex.
	 */
	private boolean hasVertex(String name) {
		for(String vertex : this.weighted.vertexSet()) {
			if (vertex.equals(name)) {
				return true;
			}
		}
		return false;
	}

       /**
         * Processes the weight part (if it exists).#
         *
         * @param tokens
         */
        void processWeight(Vector<String> tokens, String parent, String myNode) {
            if (tokens.size() == 0) {
                return;
            }
            if (tokens.get(0).equals(":")) {
                    tokens.remove(0);
                    this.weighted.setEdgeWeight(
                                    this.weighted.getEdge(parent, myNode),
                                    Double.parseDouble(tokens.get(0)));
                    tokens.remove(0);
            }
        }

	/**
	 * Parses a Newick tree.
	 * 
	 * Alters this.weighted!
	 * 
	 * The tree is passed as a set of tokens.
	 * 
	 * If some tokens are not processed, these are maintained in the vector.
	 * All consumed ones are removed.
	 * 
	 * @param tokens Stream of tokens to be parser
	 */
	void parseTree(Vector<String> tokens, String parent) throws ParseException {
		String myNode;
		if (parent == null) {
			pValue = 0;
			uuids = new Vector<String>();
		}
		//System.out.println("Top: " + parent + " ");
		//for(String tok: tokens) {
		//	System.out.print(" " + tok);
		//}
		//System.out.println();

		if (tokens.get(0).equals("(")) {
			tokens.remove(0);
			myNode = UUID.randomUUID().toString();
			uuids.add(myNode);
			this.weighted.addVertex(myNode);
			if (parent != null) {
				this.weighted.addEdge(parent, myNode);
			}
			while (!tokens.get(0).equals(")")) {
				parseTree(tokens, myNode);
				if (tokens.get(0).equals(",")) {
					tokens.remove(0);
				}
				else if (!tokens.get(0).equals(")") ) {
					throw new ParseException ("Expecting ), got " + tokens.get(0));
				}
			}
			tokens.remove(0);
			if (tokens.size() > 0) {
				String nextToken = tokens.get(0);
				char firstChar = nextToken.charAt(0);
				if (Character.isLetter(firstChar)) {
					uuids.remove(myNode);
					renameVertex(myNode, nextToken);
					myNode = nextToken;
					tokens.remove(0);
				}
				processWeight(tokens, parent, myNode);
			}
			if (parent == null) {
				String finalName;
				//Lets solve clashes
				for (String uuid: uuids) {
					do {
						finalName = this.nodePrefix + (this.pValue++);
					} while (hasVertex(finalName));
					if (uuid.equals(myNode)) {
						this.topNode = finalName;
					}
					renameVertex(uuid, finalName);
				}

				uuids = null; //for gc
			}

		}
		else {
			myNode = tokens.get(0);
			tokens.remove(0);
			this.weighted.addVertex(myNode);
			this.weighted.addEdge(parent, myNode);
			if (tokens.get(0).equals(":")) {
				tokens.remove(0);
				this.weighted.setEdgeWeight(
						this.weighted.getEdge(parent, myNode),
						Double.parseDouble(tokens.get(0)));
				tokens.remove(0);
			}
		}

	}

	/**
	 * Get given (NewickString) tree by label, converts it to weighted graph (JGraphT).
	 * 
	 * Alters this.weighted!
	 *
	 * @param label
	 * 		 label for tree selection 
	 *
	 * @return converted tree as undirectedGraph
	 */
	public WeightedGraph<String, DefaultWeightedEdge> getTreeAsWeightedJGraphT(final String label)
	throws ParseException {
		String temp;
		TreesBlock.NewickTreeString t = new TreesBlock.NewickTreeString();
	
		this.weighted =  new SimpleWeightedGraph<String, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		t = (TreesBlock.NewickTreeString) this.trees.get(label);
		temp = t.getTreeString();
		parseTree(tokenize(temp), null);

		return this.weighted;
	}
      
	/**
	 * Adds a comment.
	 * 
	 * @param comment
	 *            the comment to add.
	 */
	public void addComment(final NexusComment comment) {
		this.comments.add(comment);
	}

	/**
	 * Removes a comment.
	 * 
	 * @param comment
	 *            the comment to remove.
	 */
	public void removeComment(final NexusComment comment) {
		this.comments.remove(comment);
	}

	/**
	 * Returns all comments.
	 * 
	 * @return all the selected comments.
	 */
	public List getComments() {
		return this.comments;
	}

	protected void writeBlockContents(Writer writer) throws IOException {
		for (final Iterator i = this.comments.iterator(); i.hasNext();) {
			((NexusComment) i.next()).writeObject(writer);
			writer.write(NexusFileFormat.NEW_LINE);
		}
		writer.write(" TRANSLATE" + NexusFileFormat.NEW_LINE);
		for (final Iterator i = this.translations.entrySet().iterator(); i
				.hasNext();) {
			final Map.Entry entry = (Map.Entry) i.next();
			writer.write('\t');
			this.writeToken(writer, "" + entry.getKey());
			writer.write('\t');
			this.writeToken(writer, "" + entry.getValue());
			if (i.hasNext())
				writer.write(',');
			else
				writer.write(';');
			writer.write(NexusFileFormat.NEW_LINE);
		}
		for (final Iterator i = this.trees.entrySet().iterator(); i.hasNext();) {
			final Map.Entry entry = (Map.Entry) i.next();
			final NewickTreeString treeStr = (NewickTreeString) entry
					.getValue();
			writer.write(" TREE ");
			if (treeStr.isStarred())
				writer.write("* ");
			this.writeToken(writer, "" + entry.getKey());
			writer.write('=');
			if (treeStr.getRootType() != null)
				writer.write("[" + treeStr.getRootType() + "]");
			this.writeToken(writer, treeStr.getTreeString());
			writer.write(";" + NexusFileFormat.NEW_LINE);
		}
	}

	/**
	 * Returns the top node of the previously requested graph.
	 *
	 * The topNode will be the root node if the tree is rooted, if not
	 * then it will just be the top most node of the newick string with
	 * no biological meaning.
	 *
	 * The top node from the {@link #getTreeAsJGraphT(java.lang.String)}
	 * and {@link #getTreeAsWeightedJGraphT(java.lang.String)} might vary,
	 * and this function will return the top node of the previously called
	 * method only. If no method was called, null is returned.
	 *
	 * Note: the top node between graphs, probably does not vary, but,
	 * just in case, the note is here and the user should get a different
	 * top for each type of graph.
	 *
	 * @return the top node.
	 */
	public String getTopNode() {
		return topNode;
	}

	/**
	 * Sets the node prefix of intermediate nodes for returned graphs.
	 *
	 * The intermediate nodes of graphs have to be named, the conventions
	 * being:
	 *    1. Generate a new node name using a prefix plus an integer
	 *    2. If the node name clashes with any taxon name, use another integer
	 *
	 * @param prefix The prefix
	 */
	public void setNodePrefix(String prefix) {
		this.nodePrefix = prefix;
	}

	/**
	 * Returns the node prefix.
	 * 
	 * @return the node prefix
	 */
	public String getNodePrefix() {
		return this.nodePrefix;
	}

}
