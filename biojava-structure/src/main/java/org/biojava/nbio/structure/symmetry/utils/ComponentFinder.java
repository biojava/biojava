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
package org.biojava.nbio.structure.symmetry.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 *
 * @author Peter
 */
public class ComponentFinder<V> {
    
    private Graph<V> graph;
    private boolean[] encounter = new boolean[128];
    private List<List<V>> components;
    
    /** Creates a new instance of ComponentFinder */
    public ComponentFinder() {
    }
    
    public void setGraph(Graph<V> graph) {
        this.graph = graph;
        ensureCapacity();
        clearComponents();
        findComponents();
    }
    
    public List<V> getComponent(int index) {
        return components.get(index);
    }
    
    public List<List<V>> getComponents() {
        return components;
    }
    
    public List<V> getLargestComponent() {
        if (components.size() == 0) {
            components.add(new ArrayList<V>(0));
//            throw new NullPointerException("No components in ComponentFinder");
        }
        int largest = 0;
//        System.out.println("# components: " + components.size());
        int index = 0;
        for (int i = 0; i < components.size(); i++) {
            int size = components.get(i).size();
//            System.out.println("size: " + size);
            if (size > largest) {
                largest = size;
                index = i;
            }
        }
        return components.get(index);
    }
    
    public int getComponentCount() {
        if (components == null) return 0;
        return components.size();
    }
    
    private void findComponents() {
        // initially all vertices are marked as not encountered
        Arrays.fill(encounter, 0, graph.size(), false);
        
        for (int i = 0; i < graph.size(); i++) {
            if ( !encounter[i] ) {
                List<V> fragment = new ArrayList<V>();
                components.add(fragment);
                fragment.add(graph.getVertex(i));
                encounter[i] = true;
                traverse(i, fragment);
            }
        }
    }
    
    private void traverse(int start, List<V> fragment) {
        // mark start atom and add to fragment


        // mark all neighbors and add to fragment (recursively)
        List<Integer> neighbors = graph.getNeighborIndices(start);
        
        // if there are no neighbors, return
        if (neighbors.size() == 0) return;
        
        for (int neighbor: neighbors) {
            if (! encounter[neighbor]) {
                fragment.add(graph.getVertex(neighbor));
                encounter[neighbor] = true;
                traverse(neighbor, fragment);
            }
        }
    }
    
    private void clearComponents() {
        components = null;
        components = new ArrayList<List<V>>();
        Arrays.fill(encounter, false);
    }
    
    private void ensureCapacity() {
        if (encounter.length < graph.size()) {
            encounter = null;
            encounter = new boolean[graph.size()];
        }
    }
    
    public static void main(String[] args) {
        String va = "A";
        String vb = "B";
        String vc = "C";
        String vd = "D";
        String ve = "E";
        
//           /C-E
//        A-B |
//           \D
        SimpleGraph<String> g1 = new SimpleGraph<String>();
        g1.addVertex(va);
        g1.addVertex(vb);
        g1.addVertex(vc);
        g1.addVertex(vd);
        g1.addVertex(ve);
        g1.addEdge(va, vb);
        g1.addEdge(vb, vc);
        g1.addEdge(vc, vd);
        g1.addEdge(vd, vb);
        g1.addEdge(vc, ve);
        
        System.out.println(g1);
        
//        Fragmenter frag = new Fragmenter(g1);
//        
//        ComponentFinder c = new ComponentFinder();
////        c.setGraph(g1);
////        System.out.println("# components: " + c.getComponentCount());
//        Graph f;
//        
//        while (frag.hasNext()) {
//            f = frag.next();
//            System.out.println("frag: " + f.toString() + ", size: " + f.size());
//            c.setGraph(f);
//            System.out.println("# components: " + c.getComponentCount());
//            for (int i = 0; i < c.getComponentCount(); i++) {
//                System.out.println("Component " + c.getComponent(i));
//            }
//        }
    }
}
