package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class FlatFileCacheTest {

    final String aDNA = "ATCG";
    final String aProtein = "WCTH";

    @BeforeEach
    void before(){
        FlatFileCache.clear();
    }

    File createSmallTmpFile() throws IOException{
        File f = File.createTempFile("flatFile","txt");
       writeToFile( aDNA, f);
        return f;
    }

    private void writeToFile(String aDNA, File f) throws IOException {
        FileOutputStream fos = new FileOutputStream(f);
        fos.write(aDNA.getBytes(StandardCharsets.UTF_8));
    }

    @Test
    void flatFileRetrieve () throws IOException {
        File aDNAFile = createSmallTmpFile();
        assertEquals(0, FlatFileCache.size());
        FlatFileCache.addToCache("key", aDNAFile);
        assertEquals(1, FlatFileCache.size());

        InputStream is = FlatFileCache.getInputStream("key");
        assertNotNull(is);
        byte [] b = new byte[1024];
        int read = is.read(b);
        assertEquals(aDNAFile.length(), (long)read );
        assertEquals(aDNA, new String(b, "UTF8").substring(0,4));
    }

    @Test
    void clearRemovesAllItems () throws IOException {
        for (int i = 0; i< 10; i++) {
            FlatFileCache.addToCache(""+i, createSmallTmpFile());
        }
        assertEquals(10, FlatFileCache.size());
        FlatFileCache.clear();
        assertEquals(0, FlatFileCache.size());
    }

    @Test
    void nullReturnedIfNoValueForKey () throws IOException {
        assertNull(FlatFileCache.getInputStream("nonexistent"));
    }

    @Test
    void fileCanBeModifiedButCachedValueIsUnchanged() throws IOException{
        File aDNAFile = createSmallTmpFile();
        FlatFileCache.addToCache("key", aDNAFile);
        long originalLength = aDNAFile.length();
        
        // write new content to original file
        writeToFile( aProtein , aDNAFile);

        // retrieve from cache, is unchanged
        InputStream is = FlatFileCache.getInputStream("key");
        byte [] b = new byte[1024];
        int read = is.read(b);
        assertEquals(originalLength, (long)read );
        assertEquals(aDNA, new String(b, "UTF8").substring(0,4));
    }

}
