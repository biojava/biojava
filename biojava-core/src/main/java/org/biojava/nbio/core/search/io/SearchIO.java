/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.core.search.io;

import java.io.File;
import java.net.URL;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author pavanpa
 */
public class SearchIO implements Iterable<Result>{
    private HashMap<String,ResultFactory> extensionFactoryAssociation;
    
    /**
     * this threshold applies retrieving hsp. Those having e-value below this
     * will not be loaded.
     */
    private double evalueThreshold = Double.MAX_VALUE;
    /**
     * contains one object per query sequence describing search results.
     * Sometime also referred as Iterations.
     */
    private List<Result> results;
    
    private final String NOT_SUPPORTED_FILE_EXCEPTION = "This extension is not associated with any parser.";
    
    public SearchIO (File f)  throws Exception{
        //this(f,getFactory(f));
    }
    
    public SearchIO (File f, ResultFactory factory) throws Exception{
        factory.setFile(f);
        results = factory.createObjects(evalueThreshold);
    }
    
    public SearchIO(File f, ResultFactory factory, double maxEvalue) throws Exception{
        this.evalueThreshold = maxEvalue;
        factory.setFile(f);
        results = factory.createObjects(evalueThreshold);
    }
    
    private ResultFactory getFactory(File f){
        List<Class<?>> classes = FactoryLoader.find("Bio.SearchIO.ConcreteFactories");
        
        String filename = f.getAbsolutePath();
        int extensionPos = filename.lastIndexOf(".");
        /*
        int lastSeparator = indexOfLastSeparator(filename);
        return (lastSeparator > extensionPos ? -1 : extensionPos);
        */
        String extension = filename.substring(extensionPos + 1);
        if (extensionFactoryAssociation.get(extension) == null) 
            throw new UnsupportedOperationException(NOT_SUPPORTED_FILE_EXCEPTION
                    + "\nExtension:"+ extension);
        
        return extensionFactoryAssociation.get(extension);
    }

    public double getEvalueThreshold() {
        return evalueThreshold;
    }
    
    @Override
    public Iterator<Result> iterator() {
        return new Iterator<Result>() {
            int currentResult = 0;
            @Override
            public boolean hasNext() {
                return currentResult < results.size();
            }

            @Override
            public Result next() {
                return results.get(currentResult++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException("The remove operation is not supported by this iterator");
            }
        };
    }
    
    

    
}
/**
 * da :
 * http://stackoverflow.com/questions/15519626/how-to-get-all-classes-names-in-a-package
 * @author pavanpa
 */
class FactoryLoader {
    private static final char DOT = '.';

    private static final char SLASH = '/';

    private static final String CLASS_SUFFIX = ".class";
    
    private static final String BAD_PACKAGE_ERROR = "Unable to get resources from path '%s'. Are you sure the package '%s' exists?";
    
    public static List<Class<?>> find(String scannedPackage) {
        String scannedPath = scannedPackage.replace(DOT, SLASH);
        URL scannedUrl = Thread.currentThread().getContextClassLoader().getResource(scannedPath);
        
        if (scannedUrl == null) {
            throw new IllegalArgumentException(String.format(BAD_PACKAGE_ERROR, scannedPath, scannedPackage));
        }
        File scannedDir = new File(scannedUrl.getFile());
        List<Class<?>> classes = new ArrayList<Class<?>>();
        for (File file : scannedDir.listFiles()) {
            classes.addAll(find(file, scannedPackage));
        }
        return classes;
    }

    private static List<Class<?>> find(File file, String scannedPackage) {
        List<Class<?>> classes = new ArrayList<Class<?>>();
        String resource = scannedPackage + DOT + file.getName();
        if (file.isDirectory()) {
            for (File child : file.listFiles()) {
                classes.addAll(find(child, resource));
            }
        } else if (resource.endsWith(CLASS_SUFFIX)) {
            int endIndex = resource.length() - CLASS_SUFFIX.length();
            String className = resource.substring(0, endIndex);
            try {
                classes.add(Class.forName(className));
            } catch (ClassNotFoundException ignore) {
            }
        }
        return classes;
    }

}

