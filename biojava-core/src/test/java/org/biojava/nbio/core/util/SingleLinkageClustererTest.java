package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.Map;
import java.util.Set;

import org.junit.jupiter.api.Test;

class SingleLinkageClustererTest {
    
    // from wikipedia example
    // https://en.wikipedia.org/wiki/Single-linkage_clustering
    // it should produce clusters ((0,1),2,4),3 at distance 8.5, 10.5 and 14
    double [][] matrix = new double[][]{
        {0, 17, 21,31,23},
        {17,0,30,34,21},
        {21,30,0,28,39},
        {31,34,28,0,43},
        {23,21,39,43,0}
        };


    @Test
    void squareMatrixRequired() {
        double [][] non_square_matrix = new double[][]{{1,2},{1,2,3},{1}};
        assertThrows(IllegalArgumentException.class, ()->new SingleLinkageClusterer(non_square_matrix, false));
    }
    @Test
    void clusterWikipediaExampleDistanceMatrix(){
        SingleLinkageClusterer clusterer = new SingleLinkageClusterer(matrix, false);
        Map<Integer, Set<Integer>> result = clusterer.getClusters(Double.MAX_VALUE);
        assertEquals(5, result.get(1).size());

        result = clusterer.getClusters(0);
        assertEquals(1, result.get(1).size());
    }

    @Test
    void clusterWikipediaExampleScoreMatrix(){
        SingleLinkageClusterer clusterer = new SingleLinkageClusterer(matrix, true);
        Map<Integer, Set<Integer>> result = clusterer.getClusters(0);
        assertEquals(5, result.get(1).size());
        result = clusterer.getClusters(Double.MAX_VALUE);
        assertEquals(1, result.get(1).size());
    }

    
}
