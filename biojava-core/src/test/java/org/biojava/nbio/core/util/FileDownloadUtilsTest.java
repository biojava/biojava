package org.biojava.nbio.core.util;

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
}
