package org.biojava.nbio.core.sequence.io.embl;

import java.io.File;

public class Main {


    public static void  main(String args[]){
        File file = new File("/home/pslpt219/Desktop/embl");
        EmblParser emblParser = new EmblParser(file);
        emblParser.parse();
        System.out.println(emblParser.getSequence());
    }

}
