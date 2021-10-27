process test {

    output:
    stdout into numbers

    script:
    """
    #!/usr/bin/env python
    print("hello")
    """

}

numbers.view { "Received: $it" }