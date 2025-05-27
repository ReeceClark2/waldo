from mike.sort import TreeNode

def build_custom_tree():
    root = TreeNode("root")

    for i in range(2):  # Level 1
        l1 = TreeNode(f"Feed {i}")
        root.add_child(l1)

        for j in range(2):  # Level 2
            if j == 0:
                l2 = TreeNode(f"XX")
            else:
                l2 = TreeNode(f"YY")
            l1.add_child(l2)

            for k in range(3):  # Level 3
                if k == 0 or k == 1:
                    l3 = TreeNode(f"Calibration")
                else:
                    l3 = TreeNode(f"Data")
                l2.add_child(l3)

                # Only for first two L3s, add 2 L4 children
                if k < 2:
                    for m in range(2):
                        if m == 0:
                            l4 = TreeNode(f"Diode On")
                        else:
                            l4 = TreeNode(f"Diode Off")
                        l3.add_child(l4)
    
    return root

tree = build_custom_tree()
print(tree)
