## Modules
BioJava is composed of several submodules, one per broad bioinformatics topic. 
[biojava-core](https://github.com/biojava/biojava/tree/master/biojava-core) contains general core methods that are shared across different modules.

## Contributing
- Fork the repo
- On your fork branch off master (Ignore all other branches)
- Submit pull request from your branch to the biojava/biojava master branch
- If you are submitting a fix related to an existing issue be sure to include '#<issue number>' in the commit message to make processing your pull request easier.

## Coding
- Please use the provided formatting xml file in /development/eclipse (intellij users can import the formatter file)
- Please add javadocs following standard java conventions. Javadocs are a must for public facing API methods.
- Add `@author` tags to class javadocs.
- Be sure to add `@since` tags whenever adding a new public-facing API method/field/class.
- Removed unused imports and variables or any other obvious compiler warning.
- Please handle exceptions carefully. Follow the "Throw early, catch late" philosophy, see #111.
- Add tutorial docs to the [biojava-tutorial repository](https://github.com/biojava/biojava-tutorial)

## Testing
- We use junit 4 for testing. Some older tests are in junit 3, but all new ones should be 4.
- Unit tests should cover reasonably the code. The tests should go into the `src/test/main/` directory
of the corresponding biojava module.
- We store heavier and longer running tests into the [biojava-integrationtest module](https://github.com/biojava/biojava/tree/master/biojava-integrationtest)
- Try minimising use of external resources in tests. We haven't been very good at that in the past.
