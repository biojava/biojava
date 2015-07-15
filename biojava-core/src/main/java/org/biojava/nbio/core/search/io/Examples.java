/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.search.io;

import org.biojava.nbio.core.search.io.blast.BlastXMLQuery;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author pavanpa
 */
public class Examples {
    public static void main(String[] args) throws Exception {
        BlastXMLQuery blastXMLQuery = new BlastXMLQuery("E:\\output1668840459198813616tmp");
        LinkedHashMap<String, ArrayList<String>> hitsQueryDef = blastXMLQuery.getHitsQueryDef(1E-100);
        
        System.out.println("SearchIO test.");
        ResultFactory blastResultFactory = new BlastXMLQuery();
        SearchIO reader = new SearchIO(new File("E:\\output1668840459198813616tmp"), blastResultFactory, 1E-100);
        
        for (Result result: reader){
            System.out.println(result.getQueryDef()+"("+ result.getQueryID()+")");
            for (Hit hit: result){
                System.out.print(hit.getHitDef());
                System.out.print("(");
                for (Hsp hsp: hit){
                    System.out.print(hsp.getHspEvalue()+",");
                }
                System.out.println(")");
            }
            
        }
    }
}
