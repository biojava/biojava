package org.biojava.nbio.core.util;

import static org.biojava.nbio.core.util.FileDownloadUtils.getFileExtension;
import static org.biojava.nbio.core.util.FileDownloadUtils.getFilePrefix;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

class FileDownloadUtilsTest {

    @Nested
    class FileCopy {   
        
        private File createSrcFile () throws IOException {
            byte [] toSave = new byte []{1,2,3,4,5};
            File src = File.createTempFile("test", ".dat");
            try (FileOutputStream fos = new FileOutputStream(src);){
                fos.write(toSave);
            }
            return src;
        }

        @Test
        void copyFile() throws IOException {
            File src = createSrcFile();
            //sanity check
            assertEquals(5, src.length());
            File dest = File.createTempFile("dest", ".dat");
            assertEquals(0, dest.length());
            FileDownloadUtils.copy(src, dest);
            assertEquals(5, dest.length());

            //original is unaffected
            assertEquals(5, src.length());
            
            // bytes are identical
            try (FileInputStream fis1 = new FileInputStream(src);
                FileInputStream fis2 =  new FileInputStream(dest)) {
                int b = -1;
                while (( b = fis1.read()) != -1) {
                    int b2 = fis2.read();
                    assertEquals (b, b2);
                }
            }
        }
    }

    @Nested
    class FileExtension {
        @Test
        void getExtensionHappyCase(){
            File someFile = new File("sequence.fasta");
            assertEquals("fasta", getFileExtension(someFile));
        }

        @Test
        void lastSuffixOnlyReturned(){
            File someFile = new File("sequence.1.a.fasta");
            assertEquals("fasta", getFileExtension(someFile));
        }
        
        @Test
        void fileNameEndingInDotReturnsEmptyString(){
            File someFile = new File("noExtension.");
            assertEquals("", getFileExtension(someFile));
        }

        @Test
        void hiddenFile(){
            File someFile = new File(".m2");
            assertEquals("m2", getFileExtension(someFile));
        }

        @Test
        void noExtension(){
            File someFile = new File("nameOnly");
            assertEquals("nameOnly", getFileExtension(someFile));
        }        
    }

    @Nested
    class GetFilePrefix{
        @Test
        void standardFileName(){
            File someFile = new File("sequence.fasta");
            assertEquals("sequence", getFilePrefix(someFile));
        }
        @Test
        void prefixIsUpToFirstDot(){
            File someFile = new File("sequence.1.2.fasta");
            assertEquals("sequence", getFilePrefix(someFile));
        }

        @Test
        void noExtension(){
            File someFile = new File("nameOnly");
            assertEquals("nameOnly", getFilePrefix(someFile));
        }

        @Test
        void hiddenFile(){
            File someFile = new File(".m2");
            assertEquals("", getFilePrefix(someFile));
        }
    }

    @Nested
    class ToUnixPath {
        @Test
        void windowsToUnixAddsTrailingSlash(){
            String winPath = "C:\\a\\b\\c";
            assertEquals("C:/a/b/c/", FileDownloadUtils.toUnixPath(winPath));
        }
        @Test
        void unixPathReturnedUnchanged(){
            String path = "/a/b/c/";
            assertEquals(path, FileDownloadUtils.toUnixPath(path));
        }
    }

    @Nested
    class ToUserHome {
        String currUserHome = System.getProperty("user.home");
        @Test
        void simplePath (){
            String path="~/sequence.gb";
            assertEquals(currUserHome+File.separator+"sequence.gb", FileDownloadUtils.expandUserHome(path));
        }
        @Test
        void nestedPath (){
            String path="~/a/b/c/sequence.gb";
            assertEquals(currUserHome+File.separator+"a/b/c/sequence.gb", FileDownloadUtils.expandUserHome(path));
        }  
    }
}
