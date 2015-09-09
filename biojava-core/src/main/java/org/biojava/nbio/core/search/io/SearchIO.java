package org.biojava.nbio.core.search.io;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.ServiceLoader;

/**
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
 */

public class SearchIO implements Iterable<Result>{
    static private HashMap<String,ResultFactory> extensionFactoryAssociation;
    
    private ResultFactory factory;
    private File file;
    
    /**
     * this threshold applies in retrieving hsp. Those having e-value below this
     * will not be loaded.
     */
    private double evalueThreshold = Double.MAX_VALUE;
    /**
     * contains one object per query sequence describing search results.
     * Sometime also referred as Iterations.
     */
    private List<Result> results;
    
    private final String NOT_SUPPORTED_FILE_EXCEPTION = 
            "This extension is not associated with any parser. You can try to specify a ResultFactory object.";
    
    /**
     * Build a SearchIO reader and tries to select the appropriate parser inspecting
     * file extension.
     * 
     * @param f
     * @throws Exception 
     */
    public SearchIO (File f)  throws Exception{
        factory = guessFactory(f);
        file = f;
        process();
    }
    
    /**
     * Build a SearchIO reader and specify a ResultFactory object to be used
     * for parsing
     * 
     * @param f
     * @param factory
     * @throws Exception 
     */
    public SearchIO (File f, ResultFactory factory) throws Exception{
        file = f;
        this.factory = factory;
        process();
    }
    /**
     * Build a SearchIO reader, specify a ResultFactory object to be used for parsing
     * and filter hsp retrieved by a e-value threshold. 
     * This usually increase parsing speed
     * @param f
     * @param factory
     * @param maxEvalue
     * @throws Exception 
     */
    public SearchIO(File f, ResultFactory factory, double maxEvalue) throws Exception{
        file = f;
        this.factory = factory;
        this.evalueThreshold = maxEvalue;
        process();
    }
    
    private void process() throws Exception {
        factory.setFile(file);
        results = factory.createObjects(evalueThreshold);
    }
    
    private ResultFactory guessFactory(File f){
        if (extensionFactoryAssociation == null){
            extensionFactoryAssociation = new HashMap<String, ResultFactory>();
            ServiceLoader<ResultFactory> impl = ServiceLoader.load(ResultFactory.class);
            for (ResultFactory loadedImpl : impl) {
                List<String> fileExtensions = loadedImpl.getFileExtensions();
                for (String ext: fileExtensions) extensionFactoryAssociation.put(ext, loadedImpl);
            }
        }
               
        String filename = f.getAbsolutePath();
        int extensionPos = filename.lastIndexOf(".");
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
