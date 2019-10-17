# Releasing Pygac

1. checkout master
2. pull from repo
3. run the unittests
4. run `loghub`.  Replace <github username> and <previous version> with proper
   values.  To get the previous version run `git tag` and select the most
   recent with highest version number.

```
loghub pytroll/pygac -u <github username> -st v<previous version> -plg bug "Bugs fixed" -plg enhancement "Features added" -plg documentation "Documentation changes" -plg backwards-incompatibility "Backwards incompatible changes"
```

This command will create a CHANGELOG.temp file which need to be added
to the top of the CHANGLOG.md file.  The same content is also printed
to terminal, so that can be copy-pasted, too.  Remember to update also
the version number to the same given in step 5. Don't forget to commit
CHANGELOG.md!

5. Bump up the version with

```
bumpversion <level>
```

where level is one of `major`, `minor`, or `patch`.


6. push changes to github `git push --follow-tags`
7. Verify travis tests passed
8. Deploy sdist (and wheel) to PyPI
