package org.biojava.nbio.core.util;

import static org.biojava.nbio.core.util.FileDownloadUtils.getFileExtension;
import static org.biojava.nbio.core.util.FileDownloadUtils.getFilePrefix;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;

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
    class ExpandUserHome {
        String currUserHome = System.getProperty("user.home");
        @Test
        void minimalPath (){
        	String path="~";
        	assertEquals(currUserHome, FileDownloadUtils.expandUserHome(path));
        }
        @Test
        void simplePath (){
            String path="~/sequence.gb";
            assertEquals(currUserHome+File.separator+"sequence.gb", FileDownloadUtils.expandUserHome(path));
        }
        @Test
        void nestedPath (){
            String path="~/a/b/c/sequence.gb";
            assertEquals(currUserHome+File.separator
            		+ "a" + File.separator 
            		+ "b" + File.separator 
            		+ "c" + File.separator 
            		+ "sequence.gb", 
            		FileDownloadUtils.expandUserHome(path));
        }  
    }

    @Nested
    class URLMethods {
        final String availableUrl = "https://www.google.com";

        @Test
        void pingGoogleOK(){
            assertTrue(FileDownloadUtils.ping(availableUrl, 1000));
        }

        @Test
        void pingNonExistentFalse(){
            assertFalse(FileDownloadUtils.ping("https://non-existent.biojava", 1));
        }
    }
    @Nested
    class DeleteDirectory {

        private File createDirectoryTree () throws IOException {

            File tmpdir = Files.createTempDirectory("tmpDirPrefix").toFile();
            File child1 = new File(tmpdir, "a");
            File child2 = new File(child1, "b");
            File child3 = new File(child2, "c");
            File f = new File(child3, "seq.fa");
            child3.mkdirs();
            f.createNewFile();
            return tmpdir;
        }

        @Test
        void deleteFolderTree() throws IOException{
            File toDelete = createDirectoryTree();
            assertTrue(toDelete.exists());

            FileDownloadUtils.deleteDirectory(toDelete.getAbsolutePath());
            assertFalse(toDelete.exists());
        }
    }
}
