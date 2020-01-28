public class Configuration
 2 {
 3   private HashMap<String, String> keyValuePairs;
 4   
 5   public Configuration()   
 6   {
 7     keyValuePairs = new HashMap<String, String>();
 8     // load config from file
 9   }
10   
11   public HashMap<String, String> getValues()
12   {
13     return keyValuePairs;
14   }
15   
16   public String getValue(String key)
17   {
18     if (keyValuePairs.containsKey(key))
19       return keyValuePairs.get(key);
20     else 
21       return null;
22   }
23   
24   public void setValue(String key, String value)
25   {
26     if (keyValuePairs.containsKey(key))
27       keyValuePairs.replace(key, value);
28     else
29       keyValuePairs.put(key, value);
30     writerConfig();
31   }
32 }
